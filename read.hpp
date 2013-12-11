#include <cctype>

namespace std
{

struct read
{
public:
    // SAM fields:
    std::string QNAME;
    int FLAG;
    std::string RNAME;
    int POS;
    int MAPQ;
    std::string CIGAR;
    std::string RNEXT;
    int PNEXT;
    int TLEN;
    std::string SEQ;
    std::string QUAL;
    std::string MISC;

    // other flags:
    bool isReverse;
    int EndPOS;
    bool keep;

    int CigarM;
    int CigarD;
    int CigarI;
    int LengthOnGenome;

    int maxContigDel;

    read(std::string line)
    {
        line.append("\t");
        uint64_t pos = line.find('\t');
        uint32_t token = 0;
        keep = true;

        while (pos != std::string::npos)
        {
            ++token;

            switch (token)
            {
            case 1:
                QNAME = line.substr(0, pos);
                break;
            case 2:
                FLAG = atoi(line.substr(0, pos).c_str());
                isReverse = (0x10 & FLAG);
                break;
            case 3:
                RNAME = line.substr(0, pos);
                break;
            case 4:
                POS = atoi(line.substr(0, pos).c_str());
                break;
            case 5:
                MAPQ = atoi(line.substr(0, pos).c_str());
                break;
            case 6:
                CIGAR = line.substr(0, pos);
                break;
            case 7:
                RNEXT = line.substr(0, pos);
                break;
            case 8:
                PNEXT = atoi(line.substr(0, pos).c_str());
                break;
            case 9:
                TLEN = atoi(line.substr(0, pos).c_str());
                break;
            case 10:
                SEQ = line.substr(0, pos);
                break;
            case 11:
                QUAL = line.substr(0, pos);
                break;
            case 12:
                MISC = line.substr(0, pos);
                break;
            default:
                MISC.append("\t");
                MISC.append(line.substr(0, pos));
            }

            line.erase(0, pos+1);
            pos = line.find('\t');
        }

        calculateLengthOnGenome();
        EndPOS = POS + LengthOnGenome;
    }

    friend ostream& operator<< (ostream &out, const read& Read)
    {
        // SAM output fields
        out << Read.QNAME << '\t'
            << Read.FLAG << '\t'
            << Read.RNAME << '\t'
            << Read.POS << '\t'
            << Read.MAPQ << '\t'
            << Read.CIGAR << '\t'
            << Read.RNEXT << '\t'
            << Read.PNEXT << '\t'
            << Read.TLEN << '\t'
            << Read.SEQ << '\t'
            << Read.QUAL << '\t'
            << Read.MISC;
        return out;
    }

    void calculateLengthOnGenome()
    {
        CigarM = 0;
        CigarI = 0;
        CigarD = 0;

        char cChar;
        int cLength, start = 0;
        maxContigDel = 0;

        for (uint32_t i = 1; i < CIGAR.length(); ++i)
        {
            cChar = CIGAR.at(i);
            if (std::isalpha(cChar))
            {
                cLength = atoi(CIGAR.substr(start, i - start).c_str());
                switch (cChar)
                {
                case 'M':
                    CigarM += cLength;
                    break;
                case 'I':
                    CigarI += cLength;
                    break;
                case 'D':
                    CigarD += cLength;
                    if (cLength > maxContigDel)
                        maxContigDel = cLength;
                    break;
                }

                start = i+1;
            }
        }

        LengthOnGenome = CigarM+CigarD;
    }

    void clipNBasesFromFront(int n_bases)
    {
        int curPos, toCut, overshoot, start, cLength = 0;
        char cCigar;
        std::stringstream newCIGAR;

        curPos = this->POS;
        toCut = 0;
        start = 0;
        cCigar = '-';

        int i = 0;
        while ((curPos <= this->POS + n_bases) || (cCigar != 'M'))
        {
            cLength = 0;
            while (isdigit(cCigar = this->CIGAR.at(++i)));
            cLength = atoi(this->CIGAR.substr(start, i - start).c_str());

            if (cCigar != 'I')
            {
                curPos += cLength;
            }
            if (cCigar != 'D')
            {
                toCut += cLength;
            }

            start = i+1;
        }
        overshoot = (curPos - (this->POS + n_bases));
        int startOffset = (overshoot >= cLength ? cLength : overshoot);

        toCut -= startOffset;

        newCIGAR.clear();
        newCIGAR.str(std::string());
        newCIGAR << startOffset << cCigar << this->CIGAR.substr(start);

        /*
        // C++11:
        newCIGAR = std::to_string(startOffset);
        newCIGAR += cCigar;
        newCIGAR.append(Read.CIGAR.substr(start));
        */

        this->CIGAR = newCIGAR.str();
        this->POS = curPos - startOffset;
        this->SEQ.erase(0, toCut);
        if (this->QUAL != "*")
            this->QUAL.erase(0, toCut);
    }

    void clipNBasesFromBack(int n_bases)
    {
        std::string tempCigar;
        int curPos, toCut, overshoot, start, cLength = 0;
        char cCigar;
        std::stringstream newCIGAR;

        curPos = this->EndPOS;
        tempCigar = 'M' + this->CIGAR;
        toCut = 0;
        start = tempCigar.length() - 1;
        cCigar = tempCigar.at(start);

        int i = tempCigar.length() - 2;
        while ((curPos >= this->EndPOS - n_bases) || (cCigar != 'M'))
        {
            cCigar = tempCigar.at(start);
            cLength = 0;
            while (isdigit(tempCigar.at(--i)));
            cLength = atoi(tempCigar.substr(i+1, i+1 - start).c_str());

            if (cCigar != 'I')
            {
                curPos -= cLength;
            }
            if (cCigar != 'D')
            {
                toCut += cLength;
            }

            start = i;
            --i;
        }

        overshoot = (this->EndPOS - n_bases - curPos);
        int endOffset = (overshoot >= cLength ? cLength : overshoot);

        toCut -= endOffset;

        newCIGAR.clear();
        newCIGAR.str(std::string());
        newCIGAR << tempCigar.substr(1, start) << endOffset << cCigar;

        /*
        // C++11:
        newCIGAR = tempCigar.substr(1, start);
        newCIGAR.append(std::to_string(endOffset));
        newCIGAR += cCigar;
        */

        this->EndPOS = curPos + endOffset;
        this->CIGAR = newCIGAR.str();
        this->SEQ.erase(this->SEQ.length()-toCut, std::string::npos);
        if (this->QUAL != "*")
            this->QUAL.erase(this->QUAL.length()-toCut, std::string::npos);
    }
};

}

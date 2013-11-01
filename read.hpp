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

    // other useful flags:
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
        int pos = line.find('\t');
        int token = 0;
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

        for (int i = 1; i < CIGAR.length(); ++i)
        {
            cChar = CIGAR.at(i);
            if (std::isalpha(cChar))
            {
                cLength = atoi(CIGAR.substr(start, i - start).c_str());
                //std::cout << cChar << ": " << cLength << '\n';
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
        //std::cout << "======> " << lengthOnGenome << '\n';
        LengthOnGenome = CigarM+CigarD;
    }
};


}

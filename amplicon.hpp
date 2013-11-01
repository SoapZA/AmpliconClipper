namespace std
{

struct amplicon
{
    std::string NAME;
    int ampliconStart;
    int ampliconEnd;
    int insertStart;
    int insertEnd;
    bool isForward;
    bool isReverse;
    bool discardReads;

    amplicon(std::string line)
    {
        line.append("\t");
        int pos = line.find('\t');
        int token = 0;
        discardReads = false;

        while (pos != std::string::npos)
        {
            ++token;

            switch (token)
            {
            case 1:
                NAME = line.substr(0, pos);
                break;
            case 2:
                ampliconStart = atoi(line.substr(0, pos).c_str());
                break;
            case 3:
                ampliconEnd = atoi(line.substr(0, pos).c_str());
                break;
            case 4:
                insertStart = atoi(line.substr(0, pos).c_str());
                break;
            case 5:
                insertEnd = atoi(line.substr(0, pos).c_str());
                break;
            case 6:
                isForward = ((line.at(0) == '+') || (line.at(0) == '0'));
                isReverse = ((line.at(0) == '-') || (line.at(0) == '0'));
                break;
            case 7:
                if (line.at(0) == 'X')
                    discardReads = true;
                break;
            }

            line.erase(0, pos+1);
            pos = line.find('\t');
        }
    }

    int readOverlap(const std::read& Read)
    {
        return (std::min(this->ampliconEnd, Read.EndPOS) - std::max(this->ampliconStart, Read.POS));
    }

};

typedef std::vector< amplicon > amplicons_list;
typedef std::vector< amplicon* > amplicon_strands;
}

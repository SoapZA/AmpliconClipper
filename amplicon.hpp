
#include <utility>
#include <cassert>

namespace std
{

std::string reverse_complement(const std::string& sequence)
{
    std::string result;
    for (int i = sequence.length() - 1; i >=0; --i )
    {
        result += complement.find(sequence[i])->second;
    }
    return result;
}

uint32_t hamming_distance(const std::string& A, const std::string& B)
{
    uint32_t distance = 0;
    assert(A.length() == B.length());

    for (uint32_t i = 0; i < A.length(); ++i)
    {
        distance += (A[i] != B[i]);
    }

    return distance;
}

void find_best_overlap(const std::string& genome, const std::string& primer, int& pos, int& min_hamming_distance)
{
    min_hamming_distance = INT_MAX;
    int ham_distance;

    for (uint32_t i = 0; i < genome.length()-primer.length(); ++i)
    {
        ham_distance = hamming_distance(genome.substr(i, primer.length()), primer);

        if (ham_distance < min_hamming_distance)
        {
            pos = i;
            min_hamming_distance = ham_distance;
        }
    }
}

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
        uint64_t pos = line.find('\t');
        uint32_t token = 0;
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

    amplicon(const std::string& NAME_i, int ampliconStart_i, int ampliconEnd_i, int insertStart_i, int insertEnd_i, bool discardReads_i = false) :
        NAME(NAME_i),
        ampliconStart(ampliconStart_i),
        ampliconEnd(ampliconEnd_i),
        insertStart(insertStart_i),
        insertEnd(insertEnd_i),
        isForward(true),
        isReverse(true),
        discardReads(discardReads_i) {}

    int readOverlap(const std::read& Read)
    {
        return (std::min(this->ampliconEnd, Read.EndPOS) - std::max(this->ampliconStart, Read.POS));
    }

};

typedef std::vector< amplicon > amplicons_list;
typedef std::vector< amplicon* > amplicon_strands;

typedef std::pair< std::string, std::string > primers;

void read_amplicon_input(const std::string& input_file, std::amplicons_list& input_Amplicons, std::amplicon_strands& forward_Amplicons, std::amplicon_strands& reverse_Amplicons)
{
    std::ifstream amplicon_input_File(input_file.c_str());
    if (amplicon_input_File.is_open())
    {
        std::string line;

        getline(amplicon_input_File, line);
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        if (line[0] == '>')
        {
            // FASTA input
            if (!reference_provided)
            {
                std::cout << "Provide a reference genome with --ref in order to determine amplicon locations\n";
                exit(EXIT_FAILURE);
            }

            std::string temp = line;
            int length_to_copy;
            bool isForward;
            std::map< std::string, primers > primers_list;
            std::pair< std::map<std::string, primers>::iterator, bool> primers_list_iter;

            // load primers
            while (getline(amplicon_input_File, line))
            {
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                if (line[0] != '>')
                {
                    // primer sequence
                    temp.erase(0, 1);
                    length_to_copy = temp.length();
                    if (temp.find('_') != std::string::npos)
                    {
                        // primer pair with directionality
                        isForward = (temp[length_to_copy-1] == 'F');
                        length_to_copy -= 2;
                    }
                    temp.erase(length_to_copy);

                    primers_list_iter = primers_list.insert(std::pair< std::string, primers >(temp, primers()));
                    if (primers_list_iter.second)
                    {
                        // amplicon was created
                        (primers_list_iter.first)->second.first = line;
                    }
                    else
                    {
                        // amplicon already exists
                        if ((primers_list_iter.first)->second.second.empty())
                            (primers_list_iter.first)->second.second = line;
                        else
                        {
                            std::cout << "Too many primers for " << temp << "!\n";
                            exit(EXIT_FAILURE);
                        }
                    }

                }
                else
                {
                    temp = line;
                }
            }

            // load reference genome FASTA:
            std::string genome;
            std::ifstream genome_input_File(ref_input.c_str());
            if (genome_input_File.is_open())
            {
                getline(genome_input_File, line); // unimportant, contains the name of the reference genome

                while (getline(genome_input_File, line))
                {
                    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                    genome.append(line);
                }
            }
            else
            {
                std::cout << "The reference genome '" << ref_input << "' does not exist!\n";
                exit(EXIT_FAILURE);
            }
            genome_input_File.close();

            // determine position of primers and check for missing primer:
            for (std::map<std::string, primers>::iterator iter = primers_list.begin(); iter != primers_list.end(); ++iter)
            {
                if (iter->second.second.empty())
                {
                    std::cout << "Amplicon: " << iter->first << " is missing a second primer!\n";
                    exit(EXIT_FAILURE);
                }

                int primer1_pos, primer2_pos, primer1_length, primer2_length;
                int pos_primer, pos_primer_rc;
                int ham_distance, ham_distance_rc;

                /// First primer:
                primer1_length = iter->second.first.length();
                // search with normal primer
                find_best_overlap(genome, iter->second.first, pos_primer, ham_distance);
                // search with reverse complemented primer
                find_best_overlap(genome, reverse_complement(iter->second.first), pos_primer_rc, ham_distance_rc);

                if (ham_distance < ham_distance_rc)
                {
                    primer1_pos = pos_primer;
                }
                else
                {
                    primer1_pos = pos_primer_rc;
                }
                ++primer1_pos;

                /// Second primer:
                primer2_length = iter->second.second.length();
                // search with normal primer
                find_best_overlap(genome, iter->second.second, pos_primer, ham_distance);
                // search with reverse complemented primer
                find_best_overlap(genome, reverse_complement(iter->second.second), pos_primer_rc, ham_distance_rc);

                if (ham_distance < ham_distance_rc)
                {
                    primer2_pos = pos_primer;
                }
                else
                {
                    primer2_pos = pos_primer_rc;
                }
                ++primer2_pos;

                if (primer2_pos < primer1_pos)
                {
                    std::swap(primer1_pos, primer2_pos);
                    std::swap(primer1_length, primer2_length);
                }

                input_Amplicons.push_back(std::amplicon(iter->first, primer1_pos, primer2_pos+primer2_length, primer1_pos+primer1_length, primer2_pos, false));
            }
        }
        else
        {
            // AmpliconClipper input
            while (getline(amplicon_input_File, line))
                input_Amplicons.push_back(std::amplicon(line));
        }
    }
    else
    {
        std::cout << "Amplicon input file '" << input_file << "' does not exist!\n";
        exit(EXIT_FAILURE);
    }
    amplicon_input_File.close();

    for (uint32_t i = 0; i < input_Amplicons.size(); ++i)
    {
        std::cout << input_Amplicons[i].NAME << '\n';
        std::cout << "\tStart:         " << input_Amplicons[i].ampliconStart << '\n';
        std::cout << "\tInsert start:  " << input_Amplicons[i].insertStart << '\n';
        std::cout << "\tInsert stop:   " << input_Amplicons[i].insertEnd << '\n';
        std::cout << "\tStop:          " << input_Amplicons[i].ampliconEnd << '\n';
        std::cout << "\tDiscard reads: " << (input_Amplicons[i].discardReads ? "TRUE" : "FALSE") << '\n' << '\n';

        if (input_Amplicons[i].isForward)
            forward_Amplicons.push_back(&input_Amplicons[i]);
        if (input_Amplicons[i].isReverse)
            reverse_Amplicons.push_back(&input_Amplicons[i]);
    }
}
}

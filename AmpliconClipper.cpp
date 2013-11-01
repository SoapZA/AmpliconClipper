#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <climits>
#include <stdlib.h>
#include <cstring>

#include "read.hpp"
#include "amplicon.hpp"
#include "options.hpp"

std::vector<std::string> header_strings;
std::vector<std::read> read_container;
std::amplicons_list Amplicons;
std::amplicon_strands forwardAmplicons;
std::amplicon_strands reverseAmplicons;

void processReadWithAmplicon(const std::amplicon& Amplicon, std::read& Read)
{
    int curPos, toCut, overshoot, start, cLength;
    char cCigar;
    std::string newCIGAR;

    // trim the lower end
    if (std::min(Read.EndPOS, Amplicon.insertStart) - std::max(Read.POS, Amplicon.ampliconStart) > 0) {
        // there is overlap with the 5' PCR primer
        curPos = Read.POS;
        toCut = 0;
        start = 0;
        cCigar = '-';

        int i = 0;
        while ((curPos <= Amplicon.insertStart) || (cCigar == 'D')) {
            cLength = 0;
            while (isdigit(cCigar = Read.CIGAR.at(++i)));
            cLength = atoi(Read.CIGAR.substr(start, i - start).c_str());

            if (cCigar != 'I') {
                curPos += cLength;
            }
            if (cCigar != 'D') {
                toCut += cLength;
            }

            start = i+1;
        }
        overshoot = (curPos - Amplicon.insertStart);
        int startOffset = (overshoot >= cLength ? cLength : overshoot);

        toCut -= startOffset;

        newCIGAR = std::to_string(startOffset);
        newCIGAR += cCigar;
        newCIGAR.append(Read.CIGAR.substr(start));

        Read.CIGAR = newCIGAR;
        Read.POS = curPos - startOffset;
        Read.SEQ.erase(0, toCut);
        if (Read.QUAL != "*")
            Read.QUAL.erase(0, toCut);
    }

    // trim upper end
    std::string tempCigar;
    if (std::min(Read.EndPOS, Amplicon.ampliconEnd) - std::max(Read.POS, Amplicon.insertEnd) > 0) {
        // there is overlap with the 3' PCR primer
        curPos = Read.EndPOS;
        tempCigar = 'M' + Read.CIGAR;
        toCut = 0;
        start = tempCigar.length() - 1;

        int i = tempCigar.length() - 2;
        while ((curPos >= Amplicon.insertEnd) || (cCigar != 'M')) {
            cCigar = tempCigar.at(start);
            cLength = 0;
            while (isdigit(tempCigar.at(--i)));
            cLength = atoi(tempCigar.substr(i+1, i+1 - start).c_str());

            if (cCigar != 'I') {
                curPos -= cLength;
            }
            if (cCigar != 'D') {
                toCut += cLength;
            }

            start = i;
            --i;
        }

        overshoot = (Amplicon.insertEnd - curPos);
        int endOffset = (overshoot >= cLength ? cLength : overshoot);

        toCut -= endOffset;

        newCIGAR = tempCigar.substr(1, start);
        newCIGAR.append(std::to_string(endOffset));
        newCIGAR += cCigar;

        Read.CIGAR = newCIGAR;
        Read.SEQ.erase(Read.SEQ.length()-toCut, std::string::npos);
        if (Read.QUAL != "*")
            Read.QUAL.erase(Read.QUAL.length()-toCut, std::string::npos);
    }

    // set discard flag for unwanted amplicons
    if (Amplicon.discardReads)
        Read.keep = false;
}


int main(int argc, char** argv)
{
    parse_arguments(argc, argv);
    std::string line;

    // read SAM input file
    std::ifstream sam_input_File(sam_input.c_str());
    if (sam_input_File.is_open())
    {
        // Retrieve the header, store it in vector header_strings;
        getline(sam_input_File, line);
        while (line[0] == '@')
        {
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            header_strings.push_back(line);
            getline(sam_input_File, line);
        }

        do
        {
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            read_container.push_back(std::read(line));
            //std::cout << std::read(line) << '\n';
        }
        while (getline(sam_input_File, line));
    }
    else
    {
        std::cout << "The input file '" << sam_input << "' does not exist.\n";
        return 1;
    }
    sam_input_File.close();

    // read amplicon input file
    std::ifstream amplicon_input_File(amplicon_input.c_str());
    if (amplicon_input_File.is_open()) {
        while (getline(amplicon_input_File, line))
        {
            Amplicons.push_back(std::amplicon(line));
        }
    }
    amplicon_input_File.close();
    for (int i = 0; i < Amplicons.size(); ++i)
    {
        if (Amplicons[i].isForward)
            forwardAmplicons.push_back(&Amplicons[i]);
        if (Amplicons[i].isReverse)
            reverseAmplicons.push_back(&Amplicons[i]);
    }

    // Perform actual trimming
    std::amplicon_strands::const_iterator iter_start, iter_end, best_amplicon;
    int overlap, max_overlap = 0;
    for (int i = 0; i < read_container.size(); ++i) {
        if (read_container[i].isReverse) {
            // reverse strand read
            iter_start = reverseAmplicons.begin();
            iter_end = reverseAmplicons.end();
        }
        else {
            // forward strand read
            iter_start = forwardAmplicons.begin();
            iter_end = forwardAmplicons.end();
        }

        // choose the best amplicon based on overlap
        max_overlap = 0;
        overlap = 0;
        for (std::amplicon_strands::const_iterator j = iter_start; j != iter_end; ++j) {
            overlap = (*j)->readOverlap(read_container[i]);
            if (overlap > max_overlap) {
                max_overlap = overlap;
                best_amplicon = j;
            }
        }
        if (max_overlap > 0)
        {
            //std::cout << "Employing Amplicon: " << (*best_amplicon)->NAME << '\n';
            processReadWithAmplicon(*(*best_amplicon), read_container[i]);
            read_container[i].calculateLengthOnGenome();

            if (read_container[i].CigarI > max_insertions)
                read_container[i].keep = false;

            if (read_container[i].CigarD > max_deletions)
                read_container[i].keep = false;

            if (static_cast<double>(read_container[i].LengthOnGenome)/((*best_amplicon)->insertEnd - (*best_amplicon)->insertStart) < min_amplicon_coverage)
                read_container[i].keep = false;

            if (read_container[i].maxContigDel > max_cont_deletion)
                read_container[i].keep = false;
        }
    }

    // write final SAM output
    std::ofstream sam_output_File(sam_output.c_str());
    if (sam_output_File.is_open())
    {
        // header:
        for (int i = 0; i < header_strings.size(); ++i)
        {
            sam_output_File << header_strings[i] << '\n';
        }

        for (int i = 0; i < read_container.size(); ++i)
        {
            if (read_container[i].keep)
                sam_output_File << read_container[i] << '\n';
        }
    }
    sam_output_File.close();

    // write statistics
    if (write_statistics) {
        std::ofstream stats_output_File("stats.txt");
        stats_output_File << "M\tI\tD\n";
        for (int i = 0; i < read_container.size(); ++i)
        {
            if (read_container[i].keep)
                stats_output_File << read_container[i].CigarM << '\t' << read_container[i].CigarI << '\t' << read_container[i].CigarD << '\n';
        }

        stats_output_File.close();
    }

    return 0;
}

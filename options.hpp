#include <limits.h>
#include <map>

std::map<char, char> complement;

/* parameters: */
// input files:
std::string sam_input;
std::string sam_output;
std::string amplicon_input;
std::string ref_input;

// booleans:
bool write_statistics = false;
bool help_flag = false;
bool amplicon_input_is_fasta = false;
bool reference_provided = false;

int n_BasesFront = 0;
int n_BasesBack = 0;

int max_insertions = INT_MAX;
int max_deletions = INT_MAX;
int max_cont_deletion = INT_MAX;
double min_amplicon_coverage = 0;

void about()
{
    std::cout << "AmpliconClipper " << VERSION << "\n";
    std::cout << "\t by David Seifert 2013\n\n";
    std::cout << "Options:\n";
    std::cout << "\t -i   : input SAM file\n";
    std::cout << "\t -o   : output SAM file\n";
    std::cout << "\t -a   : input amplicon file\n";
    std::cout << "\t -S   : write statistics of reads\n";
    std::cout << "\t -f   : clip N bases from the 5' of a read\n";
    std::cout << "\t -b   : clip N bases from the 3' of a read\n";
    std::cout << "\t -I   : maximum insertions threshold\n";
    std::cout << "\t -D   : maximum deletions threshold\n";
    std::cout << "\t -C   : minimum coverage threshold as a fraction of the amplicon insert\n";
    std::cout << "\t--mD  : maximum adjacent deletions\n";
    std::cout << "\t--ref : provide the reference genome in order to locate amplicons via their primers\n";
}

#include <getopt.h>
static struct option long_options[] =
{
    {"mD", required_argument, 0, 1000},
    {"ref", required_argument, 0, 1001},
    {0, 0, 0, 0}
};

void parse_arguments(int argc, char** argv)
{
    int c, option_index = 0;

    while ((c = getopt_long (argc, argv, "i:o:a:Sf:b:hI:D:C:", long_options, &option_index)) != -1)
    {
        switch (c)
        {
            // SAM input:
        case 'i':
            sam_input = optarg;
            break;

            // SAM output:
        case 'o':
            sam_output = optarg;
            break;

            // Amplicon input
        case 'a':
            amplicon_input = optarg;
            break;

        case 'S':
            write_statistics = true;
            break;

        case 'f':
            n_BasesFront = atoi(optarg);
            break;

        case 'b':
            n_BasesBack = atoi(optarg);
            break;

        case 'h':
            help_flag = true;
            break;

        case 'I':
            max_insertions = atoi(optarg);
            break;

        case 'D':
            max_deletions = atoi(optarg);
            break;

        case 'C':
            min_amplicon_coverage = atof(optarg);
            break;

        case 1000:
            max_cont_deletion = atoi(optarg);
            break;

        case 1001:
            reference_provided = true;
            ref_input = optarg;
            break;

        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
            break;

        default:
            abort ();
        }
    }

    if ((help_flag) || (argc == 1))
    {
        about();
        exit(EXIT_SUCCESS);
    }

    complement['A'] = 'T';
    complement['T'] = 'A';
    complement['C'] = 'G';
    complement['G'] = 'C';
    complement['K'] = 'M';
    complement['M'] = 'K';
    complement['R'] = 'Y';
    complement['Y'] = 'R';
    complement['S'] = 'S';
    complement['W'] = 'W';
    complement['B'] = 'V';
    complement['V'] = 'B';
    complement['H'] = 'D';
    complement['D'] = 'H';
    complement['N'] = 'N';
}

#include <limits.h>
/* parameters: */
// input files:
std::string sam_input;
std::string sam_output;
std::string amplicon_input;

// booleans:
bool write_statistics = false;
bool help_flag = false;
double min_amplicon_coverage = 0;

int max_insertions = INT_MAX;
int max_deletions = INT_MAX;
int max_cont_deletion = INT_MAX;

void about()
{
    std::cout << "AmpliconClipper\n";
    std::cout << "\t by David Seifert 2013\n\n";
    std::cout << "Options:\n";
    std::cout << "\t-i : input SAM file\n";
    std::cout << "\t-o : output SAM file\n";
    std::cout << "\t-a : input amplicon file\n";
    std::cout << "\t-S : write statistics of reads\n";
    std::cout << "\t-I : maximum insertions threshold\n";
    std::cout << "\t-D : maximum deletions threshold\n";
    std::cout << "\t-C : minimum coverage threshold as a fraction of the amplicon insert\n";
    std::cout << "\t--mD : maximum adjacent deletions\n";
}

#include <getopt.h>
static struct option long_options[] =
{
    {"mD", required_argument, 0, 1000},
    {0, 0, 0, 0}
};

void parse_arguments(int argc, char** argv)
{
    int c, option_index = 0;

    while ((c = getopt_long (argc, argv, "i:o:a:ShI:D:C:", long_options, &option_index)) != -1)
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
}

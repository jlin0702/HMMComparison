#include "hmmComparison.h"

static void invalidArgs()
{
    std::cerr << "Please run the program with two HMM files: ./prog -query [query_file] -subject [subject_file]." << std::endl;
    exit(1);
}

static void usageAndDie(){
	std::cout << "Usage: ./prog [-h] [-query <query_file>]\n"
	<< "\t[-subject <subject_file>] [-s]\n"
    << "Details:\n"
    << "-h\n\tOpen usage manual\n"
    << "-query <query_file>\n\tFirst hmm file\n"
    << "-subject <subject_file\n\tSecond hmm file\n"
    << "-s\n\tOutput in format of score,evalue\n"
	;
	exit(0);
}

int main(int argc, char* argv[])
{
    if (argc <= 1)
        invalidArgs();

    const char* queryFile = NULL;
    const char* subjectFile = NULL;
    bool simple = false;
    if (strcmp(argv[1],"-h") == 0)
    {
        usageAndDie();
    }
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-') {
            std::string arg = argv[i];
            if (arg.substr(1) == "query")
            {
                i++;
                queryFile = argv[i];
            }
            else if (arg.substr(1) == "subject")
            {
                i++;
                subjectFile = argv[i];
            }
            else if (arg.substr(1) == "o")
            {
                i++;
                *stdout = *fopen(argv[i], "w");
            }
            else if (arg.substr(1) == "s")
            {
                simple = true;
            }
            else
            {
                std::cout << "Unrecognized argument: " << argv[i] << '\n';
                invalidArgs();
            }
        }
    }
    if (queryFile == NULL || subjectFile == NULL)
        invalidArgs();
    HMMComparison hmmComp(queryFile, subjectFile, simple);
    hmmComp.compare();
    return 0;
}
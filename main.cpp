#include "hmmComparison.h"

static void invalidArgs()
{
    std::cerr << "Please run the program with two HMM files: ./prog -query [query_file] -subject [subject_file]." << std::endl;
    exit(1);
}

int main(int argc, char* argv[])
{
    if (argc <= 1)
        invalidArgs();

   const char* queryFile = NULL;
   const char* subjectFile = NULL;
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
            else
            {
                std::cout << "Unrecognized argument: " << argv[i] << '\n';
                invalidArgs();
            }
        }
    }
    if (queryFile == NULL || subjectFile == NULL)
        invalidArgs();
    HMMComparison hmmComp(queryFile, subjectFile);
    hmmComp.compare();
    return 0;
}
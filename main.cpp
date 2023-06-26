#include "hmmComparison.h"

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cout << "Please run the program with two HMM files: ./prog [query_file] [subject_file]." << std::endl;
        return 1;
    }
    HMMComparison hmmComp(argv[1], argv[2]);
    hmmComp.compare();
    return 0;
}
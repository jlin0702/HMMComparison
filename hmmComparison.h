#ifndef _HMM_COMPARISON_H_
#define _HMM_COMPARISON_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <limits>

const char AMINOACIDS[20] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
const char DNA[4] = {'A', 'C', 'G', 'T'};

class HMMComparison
{
    private:

    struct Node
    {
        double score_;
        std::string type_;
        Node* backtrackPtr_;
        int i_;
        int j_;

        Node()
        {
            type_ = "";
            score_ = -1 * std::numeric_limits<float>::infinity();
            backtrackPtr_ = nullptr;
        }
    };

    struct HMM
    {
        std::string filePath_;
        std::string hmmType_;
        // 20: amino acid
        // 4: DNA
        int maxSymbols_;
        std::vector<std::vector<double>> stateTransitions_;
        std::vector<std::vector<double>> matchEmissions_;
        std::vector<std::vector<double>> insertEmissions_;
        double nullModel_[20];
        int length_;
    };

    HMM hmm1_;
    HMM hmm2_;
    double largestMMScore_;
    //finds where the largest indices are in S_mm matrix
    int largestI_;
    int largestJ_;
    //saa calculation
    std::vector<std::vector<double>> saaMM_;
    std::vector<std::vector<Node*>> mm_;
    std::vector<std::vector<Node*>> mi_;
    std::vector<std::vector<Node*>> im_;
    std::vector<std::vector<Node*>> dg_;
    std::vector<std::vector<Node*>> gd_;
    //amino acids with largest probability on each line for each hmm
    std::vector<char> largestSymbols1_;
    std::vector<char> largestSymbols2_;
    bool simple_;

    //calculates S_aa
    void calcSaa();

    //initializes matrices
    void initialize();

    //fills matrices according to research paper
    void fillTables();

    //recursive method to determine the path from largest in mm matrix
    //param:
    //  vector: store node sequence
    //  node: current node through
    void showTracking(std::vector<Node*>& v, Node* n);

    //prints important information
    void show();

    //process the file
    //param:
    //  hmm: which hmm to perform parsing result
    //  vector: corresponding vector to store largest probability amino acids
    void parseFile(HMM& hmm, std::vector<char>& largestAminoAcids);

    public:
    HMMComparison(std::string filePath1, std::string filePath2, bool simple);
    ~HMMComparison();

    //public function to run the code
    void compare();
};

#endif // __HMM_COMPARISON_H__
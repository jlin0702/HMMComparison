#include "hmmComparison.h"

//takes in a line of string
//returns a vector that separates the string delimited by space
std::vector<std::string> split(std::string s)
{
    std::vector<std::string> tokens;
    std::string token = "";
    s += " ";
    for (int i = 0; i < s.length(); ++i)
    {
        if (s[i] != ' ')
        {
            token += s[i];
        }
        else
        {
            if (token != "")
                tokens.push_back(token);
            token = "";
        }
    }    
    return tokens;
}

HMMComparison::HMMComparison(std::string filePath1, std::string filePath2, bool simple)
{
    largestMMScore_ = 0;
    largestI_ = 0;
    largestJ_ = 0;
    hmm1_.filePath_ = filePath1;
    hmm2_.filePath_ = filePath2;
    simple_ = simple;
    parseFile(hmm1_, largestSymbols1_);
    parseFile(hmm2_, largestSymbols2_);
    if (hmm1_.hmmType_ != hmm2_.hmmType_ )
    {
        std::cout << "Incompatible HMMs. Exiting .... \n";
        exit(1);
    }
    initialize();
}

HMMComparison::~HMMComparison()
{
    for (int i = 0; i < mm_.size(); ++i)
    {
        for (int j = 0; j < mm_[0].size(); ++j)
        {
            delete mm_[i][j];
            delete mi_[i][j];
            delete im_[i][j];
            delete dg_[i][j];
            delete gd_[i][j];
        }
    }
}

void HMMComparison::compare()
{
    calcSaa();
    fillTables();
    show();
    // for (int i = 0; i < mm_.size(); i++)
    // {
    //     for (int j = 0; j < mm_[0].size(); j++)
    //     {
    //         std::cout << mm_[i][j]->score_ << "\t";
    //     }
    //     std::cout << '\n';
    // }
    // for (int i = 0; i < saaMM_.size(); i++)
    // {
    //     for (int j = 0; j < saaMM_[0].size(); j++)
    //     {
    //         std::cout << saaMM_[i][j] << "\t";
    //     }
    //     std::cout << '\n';
    // }
}

void HMMComparison::show()
{
    std::vector<Node*> trackingVector;
    showTracking(trackingVector, mm_[largestI_][largestJ_]);
    //first: states
    //second: amino acid/dna
    std::vector<std::pair<char, char>> hmm1_state_seq;
    std::vector<std::pair<char, char>> hmm2_state_seq;
    for(int i = 0; i < trackingVector.size() ;i++)
    {
        char state1 = trackingVector[i]->type_[0];
        char symbol1;
        if (state1 == 'd' || state1 == 'g' || state1 == 'i')
            symbol1 = '-';
        else
            symbol1 = largestSymbols1_[trackingVector[i]->i_];
        hmm1_state_seq.push_back(std::make_pair(state1, symbol1));
        
        char state2 = trackingVector[i]->type_[1];
        char symbol2;
        if (state2 == 'd' || state2 == 'g' || state2 == 'i')
            symbol2 = '-';
        else
            symbol2 = largestSymbols2_[trackingVector[i]->j_];
        hmm2_state_seq.push_back(std::make_pair(state2, symbol2));
    }

    double eVal = (hmm1_.length_ * hmm2_.length_) / pow(2, largestMMScore_);
    if (!simple_)
    {
        std::cout << hmm1_.filePath_ << '\n';
        //prints the states
        std::cout << "States: \n";
        for (int i = 0; i < hmm1_state_seq.size(); ++i)
        {
            std::cout << hmm1_state_seq[i].first << ' ';
        }
        std::cout << '\n';
        //prints the amino acid/dna sequence
        if (hmm1_.hmmType_ == "amino")
            std::cout << "Amino acid sequence: \n";
        else if (hmm1_.hmmType_ == "DNA")
            std::cout << "DNA sequence: \n";
        for (int i = 0; i < hmm1_state_seq.size(); ++i)
        {
            std::cout << hmm1_state_seq[i].second << ' ';
        }
        std::cout << "\n\n";
        
        std::cout << hmm2_.filePath_ << '\n';
        //prints the states
        std::cout << "States: \n";
        for (int i = 0; i < hmm2_state_seq.size(); ++i)
        {
            std::cout << hmm2_state_seq[i].first << ' ';
        }
        std::cout << '\n';
        //prints the amino acid/dna sequence
        if (hmm2_.hmmType_ == "amino")
            std::cout << "Amino acid sequence: \n";
        else if (hmm2_.hmmType_ == "DNA")
            std::cout << "DNA sequence: \n";
        for (int i = 0; i < hmm1_state_seq.size(); ++i)
        {
            std::cout << hmm2_state_seq[i].second << ' ';
        }
        std::cout << "\n\n";
        
        std::cout << "Score: " << largestMMScore_ << '\n';
        std::cout << "E-value: " << eVal << '\n';
    }
    else
        std::cout << largestMMScore_ << "," << eVal << '\n';
}

void HMMComparison::parseFile(HMM& hmm, std::vector<char>& largestSymbols)
{
    std::ifstream inFile;
    inFile.open(hmm.filePath_);
    std::string str;
    while(inFile >> str)
    {
        if (str == "LENG")
        {
            inFile >> str;
            hmm.length_ = std::stoi(str);
        }
        else if (str == "ALPH")
        {
            inFile >> str;
            hmm.hmmType_ = str;
            if (str == "DNA")
                hmm.maxSymbols_ = 4;
            else if (str == "amino")
                hmm.maxSymbols_ = 20;
        }
        if (str == "d->d")
            break;
    }
    std::getline(inFile, str);
    std::getline(inFile, str);
    std::vector<std::string> tokens = split(str);
    for (int i = 1; i < 21; i++) {
        hmm.nullModel_[i-1] = exp(-1 * std::stod(tokens[i]));
    }
    std::getline(inFile, str);
    std::getline(inFile, str);
    int count = 0;
    while(std::getline(inFile, str))
    {
        std::vector<std::string> tokens = split(str);
        if (tokens[0] == "//")
            break;
        std::vector<double> hmmVector;
        double smallestLog = 99;
        int largestProbI = 0;
        if (count % 3 == 0)
        {
            for (int i = 1; i < hmm.maxSymbols_+1; ++i)
            {
                hmmVector.push_back(exp(-1 * std::stod(tokens[i])));
                if (std::stod(tokens[i]) < smallestLog)
                {
                    smallestLog = std::stod(tokens[i]);
                    largestProbI = i - 1;
                }
            }
            hmm.matchEmissions_.push_back(hmmVector);
            if (hmm.hmmType_ == "DNA")
                largestSymbols.push_back(DNA[largestProbI]);
            else if (hmm.hmmType_ == "amino")
                largestSymbols.push_back(AMINOACIDS[largestProbI]);
        }
        else if (count % 3 == 1)
        {
            for (int i = 0; i < hmm.maxSymbols_; ++i)
            {
                hmmVector.push_back(exp(-1 * std::stod(tokens[i])));
            }
            hmm.insertEmissions_.push_back(hmmVector);
        }
        else if (count % 3 == 2)
        {
            for (auto x : tokens)
            {
                if (x == "*")
                    hmmVector.push_back(0);
                else
                    hmmVector.push_back(exp(-1 * std::stod(x)));
            }
            hmm.stateTransitions_.push_back(hmmVector);
        }
        ++count;
    }
    inFile.close();
}

void HMMComparison::showTracking(std::vector<Node*>& v, Node* n)
{
    if(n != nullptr && n->score_ > 0)
    {
        //recursively go to first score 
        showTracking(v,n->backtrackPtr_);
        v.push_back(n);
    }
}

void HMMComparison::fillTables()
{
    for (int i = 1; i < hmm1_.length_; ++i)
    {
        for (int j = 1; j < hmm2_.length_; ++j)
        {
            //mm
            double temp1, temp2, temp3, temp4, temp5;
            temp1 = mm_[i-1][j-1]->score_ + std::log(hmm1_.stateTransitions_[i-1][0] * hmm2_.stateTransitions_[j-1][0]);
            temp2 = mi_[i-1][j-1]->score_ + std::log(hmm1_.stateTransitions_[i-1][0] * hmm2_.stateTransitions_[j-1][3]);
            temp3 = im_[i-1][j-1]->score_ + std::log(hmm1_.stateTransitions_[i-1][3] * hmm2_.stateTransitions_[j-1][0]);
            temp4 = dg_[i-1][j-1]->score_ + std::log(hmm1_.stateTransitions_[i-1][5] * hmm2_.stateTransitions_[j-1][0]);
            temp5 = gd_[i-1][j-1]->score_ + std::log(hmm1_.stateTransitions_[i-1][0] * hmm2_.stateTransitions_[j-1][5]);
            double largest = std::max(temp1, temp2);
            largest = std::max(largest, temp3);
            largest = std::max(largest, temp4);
            largest = std::max(largest, temp5);
            //calculates largest value with low of 0
            double largestWith0 = std::max(largest, 0.0);
            mm_[i][j]->score_ = saaMM_[i][j] + largestWith0;
            if (mm_[i][j]->score_ >= largestMMScore_)
            {
                largestMMScore_ = mm_[i][j]->score_;
                largestI_ = i;
                largestJ_ = j;
            }
            if (largest == temp1)
                mm_[i][j]->backtrackPtr_ = mm_[i-1][j-1];
            else if (largest == temp2)
                mm_[i][j]->backtrackPtr_ = mi_[i-1][j-1];
            else if (largest == temp3)
                mm_[i][j]->backtrackPtr_ = im_[i-1][j-1];
            else if (largest == temp4)
                mm_[i][j]->backtrackPtr_ = dg_[i-1][j-1];
            else if (largest == temp5)
                mm_[i][j]->backtrackPtr_ = gd_[i-1][j-1];
            //mi
            double temp6, temp7;
            temp6 = mm_[i-1][j]->score_ + std::log(hmm1_.stateTransitions_[i-1][0] * hmm2_.stateTransitions_[j][1]);
            temp7 = mi_[i-1][j]->score_ + std::log(hmm1_.stateTransitions_[i-1][0] * hmm2_.stateTransitions_[j][4]);
            mi_[i][j]->score_ = std::max(temp6, temp7);
            if (temp6 > temp7)
                mi_[i][j]->backtrackPtr_ = mm_[i-1][j];
            else
                mi_[i][j]->backtrackPtr_ = mi_[i-1][j];
            //im
            double temp8, temp9;
            temp8 = mm_[i][j-1]->score_ + std::log(hmm1_.stateTransitions_[i][1] * hmm2_.stateTransitions_[j-1][0]);
            temp9 = im_[i][j-1]->score_ + std::log(hmm1_.stateTransitions_[i][4] * hmm2_.stateTransitions_[j-1][0]);
            im_[i][j]->score_ = std::max(temp8, temp9);
            if (temp8 > temp9)
                im_[i][j]->backtrackPtr_ = mm_[i][j-1];
            else
                im_[i][j]->backtrackPtr_ = im_[i][j-1];
            //dg
            double temp10, temp11;
            temp10 = mm_[i-1][j]->score_ + std::log(hmm1_.stateTransitions_[i-1][2]);
            temp11 = dg_[i-1][j]->score_ + std::log(hmm1_.stateTransitions_[i-1][6]);
            dg_[i][j]->score_ = std::max(temp10, temp11);
            if (temp10 > temp11)
                dg_[i][j]->backtrackPtr_ = mm_[i-1][j];
            else
                dg_[i][j]->backtrackPtr_ = dg_[i-1][j];
            //gd
            double temp12, temp13;
            temp12 = mm_[i][j-1]->score_ + std::log(hmm2_.stateTransitions_[j-1][2]);
            temp13 = gd_[i][j-1]->score_ + std::log(hmm2_.stateTransitions_[j-1][6]);
            gd_[i][j]->score_ = std::max(temp12, temp13);
            if (temp12 > temp13)
                gd_[i][j]->backtrackPtr_ = mm_[i][j-1];
            else
                gd_[i][j]->backtrackPtr_ = gd_[i][j-1];
        }
    }
}

void HMMComparison::initialize()
{
    for (int i = 0; i < hmm1_.length_; ++i)
    {
        std::vector<Node*> v1,v2,v3,v4,v5;
        mm_.push_back(v1);
        mi_.push_back(v2);
        im_.push_back(v3);
        dg_.push_back(v4);
        gd_.push_back(v5);
        for (int j = 0; j < hmm2_.length_; ++j)
        {
            // Node n1, n2, n3, n4, n5;
            mm_[i].push_back(new Node);
            mm_[i][j]->type_ = "mm";
            if (i == 0)
            {
                mm_[i][j]->score_ = 0;
            }
            if (j == 0)
            {
                mm_[i][j]->score_ = 0;
            }
            mi_[i].push_back(new Node);
            mi_[i][j]->type_ = "mi";
            im_[i].push_back(new Node);
            im_[i][j]->type_ = "im";
            dg_[i].push_back(new Node);
            dg_[i][j]->type_ = "dg";
            gd_[i].push_back(new Node);
            gd_[i][j]->type_ = "gd";
            mm_[i][j]->i_ = mi_[i][j]->i_ = im_[i][j]->i_ = dg_[i][j]->i_ = gd_[i][j]->i_ = i;
            mm_[i][j]->j_ = mi_[i][j]->j_ = im_[i][j]->j_ = dg_[i][j]->j_ = gd_[i][j]->j_ = j;
        }
    }
}

void HMMComparison::calcSaa()
{
    for (int i = 0; i < hmm1_.length_; ++i)
    {
        std::vector<double> temp1;
        saaMM_.push_back(temp1);
        for (int j = 0; j < hmm2_.length_; ++j)
        {
            double totalMM = 0;
            for (int a = 0; a < hmm1_.maxSymbols_; ++a)
            {
                totalMM += hmm1_.matchEmissions_[i][a] * hmm2_.matchEmissions_[j][a] / ((hmm1_.nullModel_[a] + hmm2_.nullModel_[a]) / 2);
            }
            double saaMMNum = std::log(totalMM);
            saaMM_[i].push_back(saaMMNum);
        }
    }
}



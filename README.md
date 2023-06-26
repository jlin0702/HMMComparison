How to run:
1. Run "make"

2. Run program with hmm files
"./prog [hmm_file1] [hmm_file2]"
Run "./prog [hmm_file] [hmm_file] >[output_file]" to store output in file

Output:
For each HMM, there are the states and amino acid sequence.
The states correspond to the alignment pair states by concatenating the first and second state.
The amino acid sequence shows the amino acid at each state from the highest probability in the hmm file.
The score is calculated according to the algorithm.
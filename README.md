How to run:
1. Run "make"

<pre>Usage: ./prog [-h] [-query <query_file>]
    [-subject <subject_file>] [-s]
Details:
-h
    Open usage manual
-query <query_file>
    First hmm file
-subject <subject_file>
    Second hmm file
-s
    Output in format of score,evalue
</pre>

Output:<br>
For each HMM, there are the states and amino acid sequence.<br>
The states correspond to the alignment pair states by concatenating the first and second state.<br>
The amino acid sequence shows the amino acid at each state from the highest probability in the hmm file.<br>
The score is calculated according to the algorithm.
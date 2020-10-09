Computational Biology / Bioinformatics 
Lab 2 - README
Group 22

Alexandra Maroco 86369
André Branco 90013
Joana Alvoeiro 89469
Pedro Nunes 89525


This program runs the Smith Waterman algorithm for local sequence alignment, and outputs all possible optimal local alignments.
Our program is dependent on the biopython library, for the BLSOUM50 matrix. To install biopython, please use the command 'pip install biopython'.
Our program is also dependent on the numpy library. To install numpy, please use the command 'pip install numpy'.

To run:
- Open the terminal in the current directory
- Insert the command 'python3 smith_waterman.py'
- You will then be prompted to insert the first sequence to align, please do so in all caps and press enter
- You will then be prompted to insert the second sequence to align, please do so in all caps and press enter
- You will then be prompted to insert the desired gap penalty, please insert a positive integer, and press enter 
(if you wish for 4 points to be deducted for opening or continuing a gap, then the gap penalty you insert should be the number 4)
- The program will the output the score matrix, the best local alignment score, and all the possible optimal local alignments.

An example of interaction with our program would be:
(to find all optimal local alignments of sequences WPIWPC and IIWPI, with a gap penalty of 4)
($ indicates the command prompt, > indicates user input)

$ python3 smith_waterman.py
Input first sequence: 
> WPIWPC
Input second sequence: 
> IIWPI
Input gap penalty: 
> 4

Score Matrix:
[[ 0  0  0  0  0  0]
 [ 0  0  0 15 11  7]
 [ 0  0  0 11 25 21]
 [ 0  5  5  7 21 30]
 [ 0  1  2 20 17 26]
 [ 0  0  0 16 30 26]
 [ 0  0  0 12 26 28]]

The best local alignment score is 30.0

Optimal partial alignments:
Alignment 0:
WPI
WPI

Alignment 1:
IWP
IWP

# Sequence Alignment

The Algorithm of Sequence alignment is one of the most significant algorithms in bioinformatics. This file provides you with some instructions about how my Python code works, including what kind of data you need to prepare, what parameters need to be given and so on.

## Needleman_Wunsch.py

This code is implemented when you want to perform a global alignment of two sequences. Needleman-Wunsch algorithm is a very classical one to accomplish this job.

Needleman_Wunsch.py can return all the results of global alignment and give you the scores based on a very simple scoring rules. Optionally, you can choose whether to plot the matrix using heatmap or not.

If you just want to see the results of my code, you only need to modify the last three lines of the script: 

```python
seq1 = "AACGTACTCAAGTCT"
seq2 = "TCGTACTCTAACGAT"
nw(seq1, seq2, match=9, mismatch=-6, gap=-2, plot=True, plot_val=True)
```

Parameters explanation:

- **seq1**: The first sequence
- **seq2**: The second sequence

You can just type the string of the sequence (both nucleotides and amino acids are accepted)  like the example above, or use "file" function to import your sequence from a specific absolutely path: 

```
seq1 = file("file_path1")
seq2 = file("file_path2")
```

For you to notice, if you want to use function `file()` to import sequences, make sure that only one sequence in each file with no other information. JUST SEQUENCE!

As for our main function "nw", parameters are:

- **seq1** and **seq2** are two sequences you have prepared.
- **match**: If two sites are matched, how many scores can be obtained. Positive number is needed.
- **mismatch**: If two sites are mismatched, how many scores can be obtained. Negative number is needed.
- **gap**: If there are insertion sites or deletion sites, how many scores can be obtained. Negative number is needed.
- **plot**: Logical. If you want to plot the score matrix, code `plot=True`. 
- **plot_val**: Logical. If you choose to plot and you want to plot the actual value of score matrix, then code`plot_val=True`. For you to notice, if you choose to plot the value, then I **extremely recommend you** to modify the parameter "fontsize" in function `matrix_visual()` to make the plot perfect.

If you have problems or questions, welcome to contact with me!

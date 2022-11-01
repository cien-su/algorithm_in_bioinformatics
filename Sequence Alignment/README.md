# Sequence Alignment

The Algorithm of Sequence alignment is one of the most significant algorithms in bioinformatics. This file provides you with some instructions about how my Python code works, including what kind of data need you to prepare, what parameters need to be given and so on.

## Needleman_Wunsch.py

This code is implemented when you want to perform a whole alignment of two sequence. Needleman-Wunsch algorithm is a very classical one to accomplish this job.

When you run the code, you only need to modify the last three lines of the script. The variables "seq1" and "seq2" are two sequences you want to perform alignment. You can just type the string of the sequence (both nucleotides and amino acids are accepted) or use "file" function to import your sequence:

```python
seq1 = file(file_path)
```

The path of your file is needed.

For our main function "nw", the strings of two sequences are needed. The other parameters:

- match: If two sites are matched, how many scores can be obtained. The default is +8.
- mismatch: If two sites are mismatched, how many scores can be obtained. The default is -5.
- empty: If empty or blank site appears, how many scores can be obtained. The default is -3.
- plot: Logical. If you want to visualize the score matrix, code `plot=True`. If you choose to plot, then I extremely recommend you to modify the parameter "fontsize" in line 154 to make the plot perfect.

For you to notice, some Chinese annotations and "print" sentence exist in the script because I'm Chinese. However, little influence for you to use the script. If you have problems or questions, welcome to contact with me!

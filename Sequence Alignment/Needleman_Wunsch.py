# %% import pacakages
import numpy as np
import matplotlib.pyplot as plt
import time

# %% import sequences using file path
def file(path: str):
    """
    import sequences through file path
    path: absolutely file path
    """
    with open(path) as file_obj:
        seq = file_obj.read()
    # return sequence
    return seq

# %% initializing
def str_to_list(a: str, b: str):
    """
    'ATTGC'→['A', 'T', 'T', 'G', 'C']
    a: sequence 1
    b: sequence 2
    """
    return list(a), list(b)

def ini_matrix(l1: list, l2: list, gap):
    """
    initialize the score matrix
    l1: list of sequence 1
    l2: list of sequence 2
    gap: insertion or deletion score
    """
    # construct the initial matrix
    n1 = len(l1)
    n2 = len(l2)
    score_matrix = np.zeros((n1+1, n2+1))
    # initialize the edge scores
    for i in range(1, n1+1):
        score_matrix[i][0] += gap * i
    for j in range(1, n2+1):
        score_matrix[0][j] += gap * j
    # return matrix
    return score_matrix

# %% scoring
def score(matrix: np.array, l1: list, l2: list, match, mismatch, gap):
    """
    compute the score matrix
    matrix: initialized matrix
    l1: list of sequence 1
    l2: list of sequence 2
    match: match score
    mismatch: mismatch score
    gap: insertion or deletion score
    """
    # scoring in loop
    for i in range(1, len(l1)+1):
        for j in range(1, len(l2)+1):
            # compute scores from three directions
            from_left = matrix[i][j - 1] + gap  # from left
            from_above = matrix[i - 1][j] + gap  # from above
            if l1[i-1] == l2[j-1]:  # from diagonal line
                from_diag = matrix[i - 1][j - 1] + match  # match
            else:
                from_diag = matrix[i - 1][j - 1] + mismatch  # mismatch
            # compare and give the fina score
            matrix[i][j] = max(from_left, from_above, from_diag)
    return matrix

# %% trace back the path
def trace_back(res: np.array, l1: list, l2: list, match, mismatch, gap):
    """
    trace back and get the index of alignment results
    res: score matrix
    l1: list of sequence 1
    l2: list of sequence 2
    match: match score
    mismatch: mismatch score
    gap: insertion or deletion score
    """
    path = []    # to store all path results
    m_stack = [(len(l1), len(l2))]  # main stack
    a_stack = []  # assistant stack
    while m_stack:  # when main stack is not empty
        # check if reach the final
        if m_stack[-1] != (0, 0):
            # If not, then check if main stack is longer than the other. If yes, update the assistant stack
            if len(m_stack) > len(a_stack):
                a_stack.append([])
                row = m_stack[-1][0]
                col = m_stack[-1][1]
                if l1[row - 1] == l2[col - 1] and res[row][col] == res[row-1][col-1]+match:
                    a_stack[-1].append((row-1, col-1))
                elif res[row][col] == res[row-1][col-1]+mismatch:
                    a_stack[-1].append((row-1, col-1))
                if res[row][col] == res[row-1][col]+gap:
                    a_stack[-1].append((row-1, col))
                if res[row][col] == res[row][col-1]+gap:
                    a_stack[-1].append((row, col-1))
            # Check if the top element of assistant stack is empty list
            elif a_stack[-1] != []:
                m_stack.append(a_stack[-1].pop())
            # The top element of assistant stack is empty list. Pop the top of main stack
            elif a_stack[-1] == []:
                a_stack.pop()
                m_stack.pop()
        else:
            # if reach the final, then store the path and pop the top of main stack
            path.append(m_stack.copy())
            m_stack.pop()
    return path

# %% transfer index to name
def trace_to_base(match: list, l1: list, l2: list): 
    """
    transfer index to corresponding name and print
    match: list of index tuple
    l1: list of sequence 1
    l2: list of sequence 2
    """
    seq_match1 = []
    seq_match2 = []
    connect = []  # connection tag
    # loop transfer
    for index in range(1, len(match)):
        if match[index][0] != match[index - 1][0] and match[index][1] != match[index-1][1]:
            seq_match1.append(l1[match[index][0] - 1])
            seq_match2.append(l2[match[index][1] - 1])
            if l1[match[index][0] - 1] == l2[match[index][1] - 1]:
                connect.append("|")
            else:
                connect.append("*")
        elif match[index][0] == match[index - 1][0]:
            seq_match1.append("-")
            seq_match2.append(l2[match[index][1] - 1])
            connect.append(" ")
        elif match[index][1] == match[index - 1][1]:
            seq_match1.append(l1[match[index][0] - 1])
            seq_match2.append("-")
            connect.append(" ")
    # transfer list to string
    seq_match1 = "".join(seq_match1)
    seq_match2 = "".join(seq_match2)
    connect = "".join(connect)
    return seq_match1, seq_match2, connect

# %% visualizing
def matrix_visual(res_matrix: np.array, l1: list, l2: list, index: list, plot_val=False):
    """
    plot the score matrix
    res_matrix: score matrix
    l1: list of sequence 1
    l2: list of sequence 2
    index: list of index tuple
    plot_val: Logical. Whether plot matrix value on the heatmap or not. The default is False
    """
    # integrate all indexes
    path = set()
    for i in index:
        path = set(i).union(path)
    # score plot
    fig, ax = plt.subplots()
    visual_matrix = plt.imshow(res_matrix)
    plt.colorbar(visual_matrix)
    # set the axis
    l1.insert(0, '0')
    l2.insert(0, '0')
    ax.set_xticks(np.arange(0, len(l1)), minor=False)
    ax.set_yticks(np.arange(0, len(l2)), minor=False)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.set_xticklabels(l1)
    ax.set_yticklabels(l2)
    plt.xlabel("Sequence 1")
    plt.ylabel("Sequence 2")
    # check if plot_val is True
    if plot_val:
        for i in range(len(l2)):
            for j in range(len(l1)):
                if (j, i) in path:
                    ax.text(i-0.35, j+0.1,
                            res_matrix[j][i], fontsize=6, color="#FF49F5")
                else:
                    ax.text(i-0.35, j+0.1, res_matrix[j][i], fontsize=6)
    # tightly layout
    plt.tight_layout()
    plt.show()

# %% main function
def nw(seq1: str, seq2: str, match, mismatch, gap, plot=False, plot_val=False):
    """
    main function
    l1: list of sequence 1
    l2: list of sequence 2
    match: match score
    mismatch: mismatch score
    gap: insertion or deletion score
    plot: Logical. Whether plot heatmap of score matrix or not. The default is False
    plot_val: Logical. Whether plot matrix value on the heatmap or not. The default is False
    """
    # start time count
    print("\nStart Alignment......")
    start = time.time()
    # Get the lists of two sequences
    l1, l2 = str_to_list(seq1, seq2)
    # initialize score matrix
    score_matrix = ini_matrix(l1, l2, gap)
    # compute the score matrix
    res_matrix= score(score_matrix, l1, l2, match, mismatch, gap)
    # trace back the alignment
    path = trace_back(res_matrix, l1, l2, match, mismatch, gap)
    # stop time count
    stop = time.time()
    print("Finish!\n********************")
    # print the final results
    print("Print final results......")
    for i in path:
        i = i[::-1]
        seq_match1, seq_match2, connect = trace_to_base(i, l1, l2)
        print(seq_match1)
        print(connect)
        print(seq_match2)
        print("--------------------")
    print("%s results, scores: %s"%(len(path), res_matrix[-1][-1]))
    print("\nTime consuming: ", stop - start, "seconds")
    if plot:
        print("\nPlot the score matrix......")
        matrix_visual(res_matrix, l1, l2, path, plot_val)
    print("\nOver!!!")

# %% Application
seq1 = "AACGTACTCAAGTCT"
seq2 = "TCGTACTCTAACGAT"
nw(seq1, seq2, match=9, mismatch=-6, gap=-2, plot=True, plot_val=True)

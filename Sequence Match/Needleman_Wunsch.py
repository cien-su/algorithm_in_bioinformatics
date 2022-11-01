# %% 导入包
import numpy as np
import matplotlib.pyplot as plt
import time

# %% 导入序列
def file(path: str):
    """
    导入文件中的序列
    path: 文件路径
    """
    with open(path) as file_obj:
        seq = file_obj.read()
    # 返回序列
    return seq

# %% 初始化
def str_to_list(a: str, b: str):
    """
    将序列字符串转换为单个字符列表
    a: 序列一
    b: 序列二
    """
    return list(a), list(b)

def ini_matrix(l1: list, l2: list):
    """
    初始化罚分矩阵并计算边缘分值
    l1: 序列一列表
    l2: 序列二列表
    """
    # 获取序列长度构建初始矩阵
    n1 = len(l1)
    n2 = len(l2)
    score_matrix = np.zeros((n1+1, n2+1))
    # 初始化边缘分
    for i in range(1, n1+1):
        score_matrix[i][0] -= 3 * i
    for j in range(1, n2+1):
        score_matrix[0][j] -= 3 * j
    # 返回矩阵
    return score_matrix

# %% 计分
def score(matrix: np.array, l1: list, l2: list, match, mismatch, empty):
    """
    计算矩阵得分
    matrix: 初始化的矩阵
    l1: 序列一列表
    l2: 序列二列表
    match: 匹配得分
    mismatch: 不匹配时得分
    empty: 空位得分
    """
    # 循环计分
    for i in range(1, len(l1)+1):
        for j in range(1, len(l2)+1):
            # 计算三类分值
            from_left = matrix[i][j - 1] + empty  # 从左到右空位
            from_above = matrix[i - 1][j] + empty  # 从上到下空位
            if l1[i-1] == l2[j-1]:  # 对角线
                from_diag = matrix[i - 1][j - 1] + match  # 匹配
            else:
                from_diag = matrix[i - 1][j - 1] + mismatch  # 不匹配
            # 比较并赋分
            matrix[i][j] = max(from_left, from_above, from_diag)
    return matrix

# %% 回溯
def trace_back(matrix: np.array, l1: list, l2: list):
    """
    回溯获得匹配结果索引
    matrix: 结果矩阵
    l1: 序列一列表
    l2: 序列二列表
    """
    # 遍历回溯
    i = len(l1)
    j = len(l2)
    match = [(i, j)]
    while i > 1 and j > 1:
        # 找到三个值中的最大值
        max_val = max(matrix[i - 1][j], matrix[i][j - 1], matrix[i - 1][j - 1])
        if max_val == matrix[i - 1][j - 1]:
            match.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif max_val == matrix[i - 1][j]:
            match.append((i - 1, j))
            i -= 1
        elif max_val == matrix[i][j - 1]:
            match.append((i, j - 1))
            j -= 1
    if i == 1:  # 行到头，只往左取
        k = j - 1
        match += [(1, j) for j in range(k, 0, -1)]
    elif j == 1:  # 列到头，只往上取
        k = i - 1
        match += [(i, 1) for i in range(k, 0, -1)]
    # 逆序,返回坐标元组列表
    match = match[::-1]
    return match

# %% 将回溯结果转化为碱基
def trace_to_base(match: list, l1: list, l2: list):
    """
    将回溯获得的索引转化为对应的碱基并打印
    match: 回溯索引元组列表
    l1: 序列一列表
    l2: 序列二列表
    """
    seq_match1 = [l1[0]]
    seq_match2 = [l2[0]]
    for index in range(1, len(match)):
        if match[index][0] != match[index - 1][0] and match[index][1] != match[index-1][1]:
            seq_match1.append(l1[match[index][0] - 1])
            seq_match2.append(l2[match[index][1] - 1])
        elif match[index][0] == match[index - 1][0]:
            seq_match1.append("-")
            seq_match2.append(l2[match[index][1] - 1])
        elif match[index][1] == match[index - 1][1]:
            seq_match1.append(l1[match[index][0] - 1])
            seq_match2.append("-")
    seq_match1 = "".join(seq_match1)
    seq_match2 = "".join(seq_match2)
    return seq_match1, seq_match2

# %% 可视化
def matrix_visual(res_matrix: np.array, l1: list, l2: list, index: list):
    """
    可视化得分矩阵
    res_matrix: 得分矩阵
    l1: 序列一列表
    l2: 序列二列表
    index: 回溯的索引列表
    """
    # 可视化得分矩阵
    fig, ax = plt.subplots()
    visual_matrix = plt.imshow(res_matrix)
    plt.colorbar(visual_matrix)
    # 修改轴
    l1.insert(0, '0')
    l2.insert(0, '0')
    ax.set_xticks(np.arange(0, len(l2)), minor=False)
    ax.set_yticks(np.arange(0, len(l1)), minor=False)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.set_xticklabels(l2)
    ax.set_yticklabels(l1)
    plt.xlabel("Sequence 2")
    plt.ylabel("Sequence 1")
    # 标记路径上的得分
    for i in index:
        ax.text(i[1] - 0.35, i[0] + 0.1, res_matrix[i[0]][i[1]], fontsize=6)
    # 紧密排布
    plt.tight_layout()
    plt.show()

# %% 主函数
def nw(seq1: str, seq2: str, match=8, mismatch=-5, empty=-3, plot=False):
    """
    主函数
    seq1: 第一条序列
    seq2: 第二条序列
    match: 匹配得分,默认为+8
    mismatch: 不匹配时得分,默认为-5
    empty: 空位得分,默认为-3
    plot: 逻辑值,是否画出得分矩阵,默认不画出
    """
    # 开始计时
    print("\n--------------------")
    print("开始匹配......")
    start = time.time()
    # 获取列表
    l1, l2 = str_to_list(seq1, seq2)
    # 初始化打分矩阵
    score_matrix = ini_matrix(l1, l2)
    # 为矩阵赋分
    res_matrix = score(score_matrix, l1, l2, match, mismatch, empty)
    # 回溯
    match_index = trace_back(res_matrix, l1, l2)
    seq_match1, seq_match2 = trace_to_base(match_index, l1, l2)
    # 停止计时
    stop = time.time()
    print("匹配结束!")
    # 打印最终结果
    print("\n--------------------")
    print("打印得分矩阵......")
    print(res_matrix)
    print("\n--------------------")
    print("打印比对结果......")
    print(seq_match1)
    print(seq_match2)
    print("\n--------------------")
    print("打印计算耗时......")
    print("耗时为: ", stop - start, "秒")
    if plot:
        print("\n--------------------")
        print("可视化得分矩阵......")
        matrix_visual(res_matrix, l1, l2, match_index)
        print("结束!")

# %% 应用
seq1 = "AACGTACTCAAGTCT"
seq2 = "TCGTACTCTAACGAT"
nw(seq1, seq2, plot=True)
import matplotlib.pyplot as plt
import numpy as np

def fill_tmp_matrix_cols(seq):
    tmp_matrix = np.zeros((len(seq), len(seq)), dtype=object)
    for i, amino in enumerate(seq):
        tmp_matrix[:, i] = amino
    return tmp_matrix

def fill_tmp_matrix_rows(seq):
    tmp_matrix = np.zeros((len(seq), len(seq)), dtype=object)
    for i, amino in enumerate(seq):
        tmp_matrix[i, :] = amino
    return tmp_matrix



if __name__ == '__main__':
    seq1_string = 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC'
    seq2_string = 'ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC'

    seq1 = list(seq1_string)
    seq2 = list(seq2_string)

    tmp_matrix1 = fill_tmp_matrix_cols(seq1)
    tmp_matrix2 = fill_tmp_matrix_rows(seq2)

    matched = tmp_matrix1 == tmp_matrix2
    print(tmp_matrix1)
    print(tmp_matrix2)
    print(matched)
    plt.matshow(matched, 'DOT Plots')
    plt.xticks(range(len(seq1)), seq1)
    plt.yticks(range(len(seq2)), seq2)
    plt.show()

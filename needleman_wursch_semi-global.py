import numpy as np

class NeedlemanWurschSemiGlobal:
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.matrix = np.zeros((len(seq2)+1, len(seq1)+1), dtype=object)

    def fill_matrix(self, match_add=1, gap_penalty=1):
        for i in range(1, self.matrix.shape[0]):
            for j in range(1, self.matrix.shape[1]):
                # Check if there is a match
                match = True if self.seq1[j-1] == self.seq2[i-1] else False
                # Compute neighbours values
                if j == self.matrix.shape[1]-1:
                    up_value = self.matrix[i-1, j]
                else:
                    up_value = self.matrix[i-1, j] - gap_penalty
                if i == self.matrix.shape[0]-1:
                    left_value = self.matrix[i, j-1]
                else:
                    left_value = self.matrix[i, j-1] - gap_penalty
                diagonal_value = self.matrix[i-1, j-1] + match_add if match else self.matrix[i-1, j-1]
                # Take the max among the computed values
                self.matrix[i, j] = np.max(np.array([up_value, left_value, diagonal_value]))

    def find_alignments(self):
        new_seq1 = ''
        new_seq2 = ''
        match_string = ''
        score = 0
        start_position = (self.matrix.shape[0]-1, self.matrix.shape[1]-1)
        end_position = (0,0)
        current_position = start_position
        while(current_position != end_position):
            # print(current_position)
            row_index = current_position[0]
            col_index = current_position[1]
            up_reduction = 0 if (col_index==0 or col_index==self.matrix.shape[1]-1) else 1      # Allow gap in last col
            left_reduction = 0 if (row_index==0 or row_index==self.matrix.shape[0]-1) else 1    # Allow gap in last row
            # Check up
            if self.matrix[row_index-1, col_index] - up_reduction == self.matrix[row_index, col_index]:
                current_position = (row_index-1, col_index)
                new_seq1 += '-'
                new_seq2 += seq2[row_index-1]
                match_string += ' '
                score -= up_reduction
                continue
            # Check left
            if self.matrix[row_index, col_index-1] - left_reduction == self.matrix[row_index, col_index]:
                current_position = (row_index, col_index-1)
                new_seq1 += seq1[col_index-1]
                new_seq2 += '-'
                match_string += ' '
                score -= left_reduction
                continue
            # Check diagonal
            match_add = 1 if self.seq1[col_index-1] == self.seq2[row_index-1] else 0
            if self.matrix[row_index-1, col_index-1] + match_add == self.matrix[row_index, col_index]:
                current_position = (row_index-1, col_index-1)
                new_seq1 += seq1[col_index-1]
                new_seq2 += seq2[row_index-1]
                match_string += '|' if match_add==1 else ' '
                score += match_add
                continue
        return new_seq1[::-1], new_seq2[::-1], match_string[::-1], score


    def print_matrix(self):
        for row in self.matrix:
            for value in row:
                if value < 0:
                    print(value, '  ', end='', sep='')
                else:
                    print(' ', value, '  ', end='', sep='')
            print('\n')

if __name__ == '__main__':
    seq1 = 'ACACTGATCG'
    seq2 = 'CGT'

    nw_semi_global = NeedlemanWurschSemiGlobal(seq1, seq2)
    nw_semi_global.fill_matrix()
    nw_semi_global.print_matrix()
    align_seq1, align_seq2, match_string, score = nw_semi_global.find_alignments()

    print('----------------------------')
    print('original sequence 1: ', seq1)
    print('original sequence 2: ', seq2)
    print('----------------------------')
    print('alignment score: ', score)
    print('alignment:')
    print(align_seq1)
    print(match_string)
    print(align_seq2)

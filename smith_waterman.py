import numpy as np
from Bio.Align import substitution_matrices

blosum50 = substitution_matrices.load('BLOSUM50')

def get_matrix_score(m, i, j):
    return m[i][j]["score"]

def find_parent(i, j, parent_pos):
    #given a row and column, and a list of indices, determines a matrix cell's parents
    #if score is non zero, match/mismatch is always first on the parents list
    parent = []

    for pos in parent_pos:
        if(pos == 0):
            parent.append([-1, -1])
        elif(pos == 1):
            parent.append([i-1, j-1])
        elif(pos == 2):
            parent.append([i-1, j])
        else:
            parent.append([i, j-1])
    
    return parent

def get_score_matrix(m, top_score, row_n, col_n):
    #produces a matrix including only the scores
    score_matrix = np.zeros([row_n, col_n], dtype=int)
    max_score_cells = []

    for i in range(row_n):
        for j in range(col_n):
            score_matrix[i][j] = m[i][j]["score"]
            if(score_matrix[i][j] == top_score):
                max_score_cells.append([i, j])
    
    return score_matrix, max_score_cells

def display_alignments(alignments):
    print("Optimal partial alignments:")
    for i in range(len(alignments)):
        print("Alignment {}:".format(i))
        for j in range(len(alignments[i])):
            print(alignments[i][j][0], end='')
        print("")
        for j in range(len(alignments[i])):
            print(alignments[i][j][1], end='')
        print("\n")

def traceback(m, seq1, seq2, cells):
    #given a matrix where smith waterman has ben run, determines the optimal alignments
    alignments = []
    
    def traceback_aux(cell, alignment):
        parent = m[cell[0]][cell[1]]["parent"][0] #in case of multiple parents, favors the diagonal (match/mismatch), as opposed to opening gaps
                                                #then favors opening a gap in the 2nd sequence, then the 1st sequence

        if(m[cell[0]][cell[1]]["score"] == 0):
            alignment.reverse()
            return [alignment]

        elif(parent[0] == cell[0] - 1 and parent[1] == cell[1] - 1):
            alignment += [[seq1[cell[0] - 1], seq2[cell[1] - 1]]]
            
        elif(parent[0] == cell[0] - 1 and parent[1] == cell[1]):
            alignment += [[seq1[cell[0] - 1], "-"]]

        else:
            alignment += [["-", seq2[cell[1] - 1]]]
        
        return traceback_aux(parent, alignment)

    for cell in cells:
        new_align = traceback_aux(cell, [])
        alignments += new_align
    
    return alignments
        

def smith_waterman():
    #runs smith waterman over a matrix characterized by user input parameters, and displays results
    A = input("Input first sequence: \n")
    B = input("Input second sequence: \n")
    gp = int(input("Input gap penalty: \n"))

    rows = len(A) + 1
    cols = len(B) + 1
    matrix = np.empty([rows, cols], dtype=object)

    for i in range(rows):
        matrix[i][0] = {"score": 0, "parent": [[-1,-1]]}

    for j in range(cols):
        matrix[0][j] = {"score": 0, "parent": [[-1,-1]]}

    max_score = -1

    for i in range(1, rows):
        for j in range(1, cols):
            scores = [0, get_matrix_score(matrix, i-1, j-1) + blosum50[(A[i-1], B[j-1])], get_matrix_score(matrix, i - 1, j) - gp, get_matrix_score(matrix, i, j - 1) - gp]
            chosen_score = max(scores)
            parent_pos = [n for n in range(4) if scores[n] == chosen_score]
            matrix[i, j] = {"score": chosen_score, "parent": find_parent(i, j, parent_pos)}
            if(chosen_score > max_score):
                max_score = chosen_score

    score_matrix, max_score_cells = get_score_matrix(matrix, max_score, rows, cols)

    print("\nScore Matrix:\n{}\n".format(score_matrix))
    print("The best local alignment score is {}\n".format(max_score))

    display_alignments(traceback(matrix, A, B, max_score_cells))

smith_waterman()

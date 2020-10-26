import numpy as np

def viterbi():

    trans = np.array([[0.6,0.4,0],[0.25,0.5,0.25],[0.25,0.25,0.5]])
    emission = np.array([[0.4,0.3,0,0.3],[0.1,0.1,0.4,0.4],[0.4,0.3,0.3,0]]) # order: A-T-C-G

    S = input("Input the sequence: \n")

    v = {1: [[1,0],[2,0]], 2: [[1,0],[2,0],[3,0]], 3: [[1,0],[2,0],[3,0]]}


viterbi()import numpy as np

def viterbi():

    trans = np.array([[0.6,0.4,0],[0.25,0.5,0.25],[0.25,0.25,0.5]])
    emission = np.array([[0.4,0.3,0,0.3],[0.1,0.1,0.4,0.4],[0.4,0.3,0.3,0]]) # order: A-T-C-G

    S = input("Input the sequence: \n")

    v = {1: [[1,0],[2,0]], 2: [[1,0],[2,0],[3,0]], 3: [[1,0],[2,0],[3,0]]}


viterbi()

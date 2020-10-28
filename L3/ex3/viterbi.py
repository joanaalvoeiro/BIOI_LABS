import numpy as np
import math

def viterbi():
    states = [1, 2, 3]
    trans = np.array([[0.6,0.4,0],[0.25,0.5,0.25],[0.25,0.25,0.5]])
    emission = {1:{"a": 0.4, "t": 0.3, "c": 0, "g": 0.3}, 2:{"a": 0.1, "t": 0.1, "c": 0.4, "g": 0.4},
                3:{"a": 0.4, "t": 0.3, "c": 0.3, "g": 0}}

    S = input("Input the sequence: \n")
    S = S.lower()

    weight_graph = []
    for i in range(len(S)):#trellis creation
        node = []
        for _ in range(len(states)):
            node.append([])
        weight_graph.append(node)
    
    for i in range(len(states)):#initialize first column
        weight_graph[0][i] = [0, emission[i+1][S[0]]/3]

    for i in range(1, len(S)):#virtebi algorithm
        for state in states:
            sources = [s for s in states]
            weights = []

            for s in range(len(sources)):
                source_w = weight_graph[i - 1][s][1]
                trans_p = trans[s][state - 1]
                emission_p = emission[state][S[i]]
                if(emission_p != 0 and trans_p != 0):
                    weights.append(math.log(emission_p) + source_w + math.log(trans_p))
                else:
                    weights.append(float('-inf'))

            max_w = max(weights)
            parent = weights.index(max_w) + 1
            weight_graph[i][state - 1] = [parent, max_w]
    
    max_total_weight = float('-inf') 
    max_state = 0
    
    for i in range(len(states)):#find the largest probability
        if(weight_graph[len(S) - 1][i][1] > max_total_weight):
            max_total_weight = weight_graph[len(S) - 1][i][1]
            max_state = i + 1
    
    i = len(S) - 1
    seq = [max_state]

    curr_state = max_state
    while (i > 0):#traceback
        parent = weight_graph[i][curr_state - 1][0]
        seq.append(parent)
        curr_state = parent
        i -= 1

    seq.reverse()

    print(seq)

viterbi()
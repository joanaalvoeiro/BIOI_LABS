import numpy as np
import math

states = [1, 2, 3]
trans = np.array([[0.6,0.4,0],[0.25,0.5,0.25],[0.25,0.25,0.5]])
emission = {1:{"a": 0.4, "t": 0.3, "c": 0, "g": 0.3}, 2:{"a": 0.1, "t": 0.1, "c": 0.4, "g": 0.4},
            3:{"a": 0.4, "t": 0.3, "c": 0.3, "g": 0}}

S = input("Input the sequence: \n")
S = S.lower()

def create_trellis(sequence, trellis_states, emissions):
    trellis = []
    for i in range(len(sequence)):
        node = []
        for _ in range(len(trellis_states)):
            node.append([])
        trellis.append(node)
    
    for i in range(len(trellis_states)):#initialize first column
        trellis[0][i] = [0, emissions[i+1][sequence[0]]/len(trellis_states)]
    
    return trellis

def viterbi(v_states, v_trans, v_emission, v_seq):
    weight_graph = create_trellis(v_seq, v_states, v_emission)

    for i in range(1, len(v_seq)):#compute probabilities
        for state in v_states:
            weights = []
            for s in range(len(v_states)):
                source_w = weight_graph[i - 1][s][1]
                trans_p = v_trans[s][state - 1]
                emission_p = v_emission[state][v_seq[i]]
                if(emission_p != 0 and trans_p != 0):
                    weights.append(math.log(emission_p) + source_w + math.log(trans_p))
                else:
                    weights.append(float('-inf'))

            max_w = max(weights)
            parent = weights.index(max_w) + 1
            weight_graph[i][state - 1] = [parent, max_w]
    
    max_total_weight = float('-inf') 
    max_state = 0
    
    for i in range(len(v_states)):#find the largest probability
        if(weight_graph[len(v_seq) - 1][i][1] > max_total_weight):
            max_total_weight = weight_graph[len(v_seq) - 1][i][1]
            max_state = i + 1
    
    i = len(v_seq) - 1
    seq = [max_state]

    curr_state = max_state
    while (i > 0):#traceback
        parent = weight_graph[i][curr_state - 1][0]
        seq.append(parent)
        curr_state = parent
        i -= 1

    seq.reverse()

    print('')
    print("Pi*:", seq)
    print('')

def forward(f_states, f_trans, f_emission, f_seq):
    probabilities = np.zeros((len(f_states), len(f_seq)))

    for i in range(len(f_states)):#initialize first column
        probabilities[i][0] = f_emission[i+1][f_seq[0]]/len(f_states)

    for i in range(1, len(f_seq)):#compute probabilities
        for state in f_states:
            for s in range(len(f_states)):
                source_p = probabilities[s][i - 1]
                trans_p = f_trans[s][state - 1]
                emission_p = f_emission[state][f_seq[i]]

                probabilities[state - 1][i] += emission_p * trans_p * source_p
    
    total_p = np.sum(probabilities[:, len(f_seq) - 1])

    print("Probability matrix for Forward algorithm (states x sequence):\n", probabilities)
    print('')
    print("Probability of given sequence: ", total_p)

viterbi(states, trans, emission, S)
forward(states, trans, emission, S)
import numpy as np
import math

states = [1, 2, 3]
trans = np.array([[0.6,0.4,0],[0.25,0.5,0.25],[0.25,0.25,0.5]])
emission = {1:{"a": 0.4, "t": 0.3, "c": 0, "g": 0.3}, 2:{"a": 0.1, "t": 0.1, "c": 0.4, "g": 0.4},
            3:{"a": 0.4, "t": 0.3, "c": 0.3, "g": 0}}

S = input("Input the sequence: \n")
S = S.lower()

index = int(input("Input the index: \n"))

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

    print('')
    print("Probability matrix for Forward algorithm (states x sequence):\n", probabilities)
    print('')
    print("Probability of given sequence: ", total_p)

    return probabilities, total_p

def backward(b_states, b_trans, b_emission, b_seq):
    probabilities = np.zeros((len(b_states), len(b_seq)))

    for i in range(len(b_states)):#initialize last column
        probabilities[i][len(b_seq) - 1] = 1 #not sure????

    for i in range(len(b_seq) - 2, -1, -1):#compute probabilities
        for state in b_states:
            for s in range(len(b_states)):
                source_p = probabilities[s][i + 1]
                trans_p = b_trans[state -1][s]
                emission_p = b_emission[s + 1][b_seq[i + 1]]

                probabilities[state - 1][i] += emission_p * trans_p * source_p
    
    total_p = 0

    for i in range(len(b_states)):
        source_p = probabilities[i][0]
        trans_p = 1 / len(b_states)
        emission_p = b_emission[i + 1][b_seq[0]]

        total_p += trans_p * emission_p * source_p

    print('')
    print("Probability matrix for Backward algorithm (states x sequence):\n", probabilities)
    print('')
    print("Probability of given sequence: ", total_p)

    return probabilities

def forward_backward(fb_states, fb_trans, fb_emission, fb_seq, fb_index):
    forward_m, px = forward(fb_states, fb_trans, fb_emission, fb_seq)
    backward_m = backward(fb_states, fb_trans, fb_emission, fb_seq)

    print('')
    posterior_probs = []
    for s in range(len(fb_states)):
        fi = forward_m[s][index - 1]
        bi = backward_m[s][index - 1]
        posterior_probs.append((fi * bi) / px)
        print("P(Pi{} = {} | {}) = {}".format(fb_index, s + 1, fb_seq, posterior_probs[s]))

    posterior_prob = max(posterior_probs)

    print('')
    print("Pi^{} = {}".format(fb_index, posterior_prob))

forward_backward(states, trans, emission, S, index)
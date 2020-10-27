import numpy as np

def viterbi():
    states = [1, 2, 3]
    trans = np.array([[0.6,0.4,0],[0.25,0.5,0.25],[0.25,0.25,0.5]])
    emission = {1:{"a": 0.4, "t": 0.3, "c": 0, "g": 0.3}, 2:{"a": 0.1, "t": 0.1, "c": 0.4, "g": 0.4},
                3:{"a": 0.4, "t": 0.3, "c": 0.3, "g": 0}}

    #S = input("Input the sequence: \n")
    S="CATGCGGGTTATAAC"
    S = S.lower()

    graph = [[] for _ in range(len(S))]

    graph[0].append([[[0, emission[1][S[0]]]], [[0, emission[2][S[0]]]], [[0, emission[3][S[0]]]]])

    for i in range(1, len(S)):
        for state_init in states:
            curr_node = []
            for state_final in states:
                weight =  emission[state_final][S[i]] * trans[state_init - 1][state_final - 1]
                curr_node.append([state_final, weight])
            graph[i].append(curr_node)  
    
    weight_graph = []
    for i in range(len(S)):
        node = []
        for _ in range(len(states)):
            node.append([])
        weight_graph.append(node)
    
    for i in range(len(states)):
        weight_graph[0][i] = [0, emission[i+1][S[0]]]
    

    for i in range(1, len(S)):
        for state in states:
            edges = [e for e in graph[i][state - 1]]
            weights = []
            for e in range(len(edges)):
                weights.append(edges[e][1])
            max_w = max(weights)
            parent = weights.index(max_w) + 1
            weight_graph[i][state - 1] = [parent, max_w]
    

    max_total_weight = 0
    max_state = 0
    
    for i in range(len(states)):
        if(weight_graph[len(S) - 1][i][1] > max_total_weight):
            max_total_weight = weight_graph[len(S) - 1][i][1]
            max_state = i + 1
    
    i = len(S) - 1
    seq = [max_state]

    curr_state = max_state
    while (i > 1):
        parent = weight_graph[i][curr_state][0]
        seq.append(parent)
        i -= 1

    seq.reverse()

    print(seq)

viterbi()
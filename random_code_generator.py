import random
import numpy as np
import itertools
from numpy.random import rand

from pyldpc.code import parity_check_matrix

def generate_random_parity_check(n,dv,dc):

    if n*dv % dc != 0:
        return [0]

    k = int(n*(1-dv/dc))


    # Generate a list containing each check node dc times which will be used to draw from
    check_nodes_available = []
    for i in range(n-k):
        check_nodes_available += dc*[i]

    # Add try/except to handle the case that check_nodes_available_to_variable becomes empty
    # in which case, random.choice throws an error
    try:
        variable_nodes = {}
        for i in range(n):
            variable_nodes[i] = []
            check_nodes_available_to_variable = check_nodes_available.copy()
            for selction in range(dv):
                choice = random.choice(check_nodes_available_to_variable)
                variable_nodes[i] += [choice]
                check_nodes_available.remove(choice)
                check_nodes_available_to_variable = [x for x in check_nodes_available_to_variable if x != choice]
    except:
        return generate_random_parity_check(n,dv,dc)
    

    # Check that all variable nodes have all dv check nodes connected 
    # This does not always happen because when you get to the last variable node there
    # may be no check nodes available without repeating one it already has
    for check_node_combinations in variable_nodes.values():
        if len(check_node_combinations) != dv:
            return generate_random_parity_check(n,dv,dc)
        check_node_combinations.sort()

    # Check that no two lists of check nodes are identical for any variable node 
    # (this would result in two identical rows in the parity check matrix)
    for check_node_combination in variable_nodes.values():
        count = 0
        for other_check_node_combinations in variable_nodes.values():
            if other_check_node_combinations == check_node_combination:
                count += 1
        if count > 1:
            return generate_random_parity_check(n,dv,dc)

    # Convert dictionary of variable nodes to parity check matrix
    parity_check_matrixT = np.zeros((n, n-k))
    for variable_node, connected_check_nodes in variable_nodes.items():
        parity_check_matrixT[variable_node][connected_check_nodes] = 1
    return np.array(parity_check_matrixT, dtype='int0').T

def generate_random_parity_check_no_checks(n, dv, dc):
    k = int(n*(1-dv/dc))

    # Randomise 'socket' connections
    random_vector = np.random.permutation(n*dv)

    # Divide each socket number by dv and round down, to get the values of each variable node
    random_vector = random_vector // dv
    random_vector = np.array(random_vector, dtype='int')

    random_vector = np.reshape(random_vector, (n-k, dc))

    parity_check_matrix = np.empty((n-k,n), dtype='bool')
    
    for index, check_node in enumerate(random_vector):
        if len(np.unique(check_node)) == dc:
            all_zero = np.zeros(n)
            np.put(all_zero, check_node, True)
            parity_check_matrix[index] = all_zero
        else:
            return generate_random_parity_check_no_checks(n,dv,dc)
    
    return parity_check_matrix

for i in range(1000):
    print(i)
    generate_random_parity_check_no_checks(10000,3,6)

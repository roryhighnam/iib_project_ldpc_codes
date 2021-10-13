import random
import numpy as np

def generate_random_parity_check(n,dv,dc):

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

# for i in range(10):
#     print(generate_random_parity_check(10,5,3,6))

print(generate_random_parity_check(500,3,6))
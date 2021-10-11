import numpy as np
import galois
from pyldpc.code import coding_matrix_systematic, coding_matrix
from pyldpc.utils import gausselimination
from channels import BEC
from random_code_generator import generate_random_parity_check

class regular_LDPC_code():

    def __init__(self,parity_check):
        self.generator_matrix = coding_matrix(parity_check)
        self.parity_check = parity_check
        self.n = len(parity_check[0])
        self.k =  len(parity_check[0]) - len(parity_check)
        self.dv = np.sum(parity_check, axis=0)[0]
        self.dc = np.sum(parity_check, axis=1)[0]
        self.rate = self.k / self.n

    def encode(self,binary_sequence):
        '''Use generator matrix to encode the input binary sequence'''
        assert len(binary_sequence) == self.k, "The input sequence must have the correct input length (equal to k)"
        return np.mod(np.dot(self.generator_matrix, binary_sequence),2)

    def optimal_decode(self, binary_sequence):
        '''Use Gaussian Elimination and back subsitution to solve for the erased bits,
        given the values of parity check equations for the non-erased bits. If solution cannot
        be found, returns -1.'''
        # remaining_parity_checks will contain the columns of the parity check matrix where there are erasures
        remaining_parity_checks = []
        # Conversely, known_parity_checks will contain the columns of the parity check matrix where there are not erasures
        known_parity_checks = []
        # known_codeword contains the received bits that are not erasures
        known_codeword = []

        # Iterate through the input sequence and add corresponding data to variables above
        for index,bit in enumerate(binary_sequence):
            if bit != 0:
                # If bit is not erasure, add column of parity check matrix to known_parity_checks for corresponding column. Also add the bit value to known_codeword.
                known_parity_checks += [self.parity_check[:,index]]
                if bit == -1:
                    known_codeword = np.append(known_codeword, 0)
                else:
                    known_codeword = np.append(known_codeword, bit)
            else:
                # If bit is erasure, add the corresponding column to remaining_parity_check
                remaining_parity_checks += [self.parity_check[:,index]]

        known_parity_checks = np.array(known_parity_checks).T
        remaining_parity_checks = np.array(remaining_parity_checks).T
        known_codeword = np.array(known_codeword)
        # target is the values of the parity check equations for the known bits (i.e. non-erasures)
        target = np.mod(np.dot(known_parity_checks, known_codeword),2)

        # Deal with underdetermined matrix or if there are no erasures
        if remaining_parity_checks.shape[0] < binary_sequence.count(0):
            # Should we return corrected bits if we can solve for some bits?
            return [-1]
        if len(known_codeword) == self.n:
            return list(map(lambda x: 0 if x==-1 else x, binary_sequence))
        
        # Perform gaussian elimination to recover decoded codeword (This could still fail if we cannot use back substitution)
        upper_triangular_matrix, target = gausselimination(remaining_parity_checks, target)
        try:
            solved_unknowns = np.linalg.solve(upper_triangular_matrix[:binary_sequence.count(0),:], target[:binary_sequence.count(0)])
        except:
            return [-1]

        # Finally, loop through input sequence and convert -1s back to 0s and add in our solved unknowns (previously erasures)
        decoded_codeword = []
        for bit in binary_sequence:
            if bit != 0:
                if bit == -1:
                    decoded_codeword += [0]
                else:
                    decoded_codeword += [1]
            else:
                decoded_codeword += [solved_unknowns[0]]
                solved_unknowns = np.delete(solved_unknowns, 0)

        return np.mod(np.array(decoded_codeword, dtype='int0'),2)

    def compute_cv_message(self, inputs):
        '''Function to compute the check->variable message. Pass in all variable->check 
        messages except from the variable you are sending message to'''
        if 0 in inputs:
            return 0
        else:
            modified_inputs = list(map(lambda x: 0 if x==-1 else x, inputs))
            if sum(modified_inputs)%2 == 0:
                return -1
            else:
                return 1
    
    def compute_vc_message(self, inputs):
        '''Function to compute v->c message. Pass in all check->variable messages
        except from the check you are sending message to.'''
        if -1 in inputs:
            return -1
        elif 1 in inputs:
            return 1
        else:
            return 0

    def message_pass_decode(self, binary_sequence, max_its):
        # Generate lookup tables for connected variable and check nodes ie dictionary 
        # with each check as key and connected variables list as value
        check_lookup = {}
        for row_index,row in enumerate(self.parity_check):
            check_list = []
            for index,value in enumerate(row):
                if value==1:
                    check_list.append(index)
            check_lookup[row_index] = check_list

        variable_lookup = {}
        for col_index,col in enumerate(self.parity_check.T):
            var_list = []
            for index,value in enumerate(col):
                if value==1:
                    var_list.append(index)
            variable_lookup[col_index] = var_list

        # Matrices that contain all message updates - not very memory efficient for larger codes
        Mcv = np.zeros_like(self.parity_check)
        Mvc = np.zeros_like(self.parity_check)

        # Initialise v->c messages as the channel output for each variable v
        for variable_node,check_node_list in variable_lookup.items():
            for check in check_node_list:
                Mvc[check][variable_node] = binary_sequence[variable_node]

        # Iterate for a fixed number of iterations, specified by max_its
        # May want to change this to stop iterating when messages do not update?
        for it in range(max_its):
            # Update Mcv first:
            for check_node, variable_node_list in check_lookup.items():
                for variable_node in variable_node_list:
                    # Remove current variable node from inputs to compute_cv_message function
                    other_variable_nodes = variable_node_list.copy()
                    other_variable_nodes.remove(variable_node)
                    Mcv[check_node][variable_node] = self.compute_cv_message(Mvc[check_node][other_variable_nodes])

            # Then update Mvc
            for variable_node, check_node_list in variable_lookup.items():
                for check_node in check_node_list:
                    # If current Mvc message is already non-erasure, do nothing
                    if Mvc[check_node][variable_node] == 0:
                        # Otherwise, run compute_vc_message with check_nodes other than the current one
                        other_check_nodes = check_node_list.copy()
                        other_check_nodes.remove(check_node)
                        Mvc[check_node][variable_node] = self.compute_vc_message(Mcv.T[variable_node][other_check_nodes])
        # Convert messages to a list of decoded values. On the BEC there are no conflicts so we can simply sum each column
        # and divide by variable degree as it is a regular code.
        decoded_values = np.array(np.sum(Mvc, axis=0)/self.dv, dtype='int0')
        # Convert decoded values back to binary list, leaving '?' if erasures still present
        decoded_values = list(map(lambda x:'?' if x==0 else x, decoded_values))
        return list(map(lambda x: 0 if x==-1 else x, decoded_values))

# Example taken from lecture notes (n=20, k=10, dv=3, dc=6)
# parity_check = np.array([[0,0,1,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,1,0],[0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0],[1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,1,0],
#                         [0,1,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1],[1,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0],[0,0,1,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1],
#                         [0,0,0,0,1,1,0,1,0,1,0,1,0,0,1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,0],[0,1,0,1,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,1],
#                         [0,0,1,0,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0]])
# test_sequence = [0,1,1,0,0,0,1,1,0,1]

# Rank of parity check can be used to test for independent rows (or equivalently redundant parity check equations)
# GF2 = galois.GF(2)
# print(np.linalg.matrix_rank(GF2(parity_check)))

# Example of bad parity check matrix - not full rank (rank of 3)
# parity_check = np.array([[0,1,0,1,1,0,0,1],[1,1,1,0,0,1,0,0],[0,0,1,0,0,1,1,1],[1,0,0,1,1,0,1,0]])
# test_sequence = [1,0,0,1]
# print(np.linalg.matrix_rank(GF2(parity_check)))

# parity_check = np.array([[1,0,1,1,1,0,1,1,0,0],[0,1,1,0,0,1,1,0,1,1],[1,1,0,1,1,1,0,1,0,0],[0,1,1,0,1,0,1,0,1,1],[1,0,0,1,0,1,0,1,1,1]])
# test_sequence = [1,1,0,1,0]
# print(np.linalg.matrix_rank(GF2(parity_check)))

parity_check = generate_random_parity_check(20,10,3,6)

LDPC = regular_LDPC_code(parity_check)

test_sequence = np.random.randint(2, size=LDPC.k)

BEC = BEC(0.3)
print('Generated codeword:', LDPC.encode(test_sequence))

channel_output = BEC.transmit(LDPC.encode(test_sequence))
print('Channel output:', channel_output)

print('Number of erasures:', channel_output.count(0))

decoded_codeword = LDPC.optimal_decode(channel_output)
print('Optimally decoded codeword:', decoded_codeword)

print('Message passing codeword:', LDPC.message_pass_decode(channel_output, 50))

for i in range(50):
    # print('Generated Codeword:', LDPC.encode(test_sequence))
    channel_output = BEC.transmit(LDPC.encode(test_sequence))
    # print('Channel output:', channel_output)
    decoded_codeword = LDPC.message_pass_decode(channel_output, 50)
    # print('Decoded Codeword:', decoded_codeword)
    if '?' in decoded_codeword:
        print('Still contains errors!')
        if -1 in LDPC.optimal_decode(channel_output):
            print('Optimal decoder successful')
        else:
            print('Optimal decoder also cannot solve')
        print()
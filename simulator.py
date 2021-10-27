from random import gauss
import numpy as np
import galois
from pyldpc.code import coding_matrix_systematic, coding_matrix
from pyldpc.utils import gausselimination, gaussjordan
from channels import BEC
from random_code_generator import generate_random_parity_check
import matplotlib.pyplot as plt
import csv
import multiprocessing
import cProfile
import re

class regular_LDPC_code():

    def __init__(self,parity_check):
        # self.generator_matrix = coding_matrix(parity_check)
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
        # Deal with no erasures case
        if len(known_codeword) == self.n:
            return list(map(lambda x: 0 if x==-1 else x, binary_sequence))
        
        # Perform gaussian elimination to recover decoded codeword
        # Check all diagonals have a 1, otherwise remove this column -> this bit cannot be determined

        upper_triangular_matrix, target = gausselimination(remaining_parity_checks, target)

        columns_to_delete = []
        for i in range(upper_triangular_matrix.shape[1]):
            if upper_triangular_matrix[i][i] != 1:
                columns_to_delete.append(i)
        upper_triangular_matrix = np.delete(upper_triangular_matrix, columns_to_delete, axis=1)
        upper_triangular_matrix = np.delete(upper_triangular_matrix, columns_to_delete, axis=0)
        target = np.delete(target, columns_to_delete)

        try:
            solved_unknowns = np.mod(np.linalg.solve(upper_triangular_matrix[:upper_triangular_matrix.shape[1],:], target[:upper_triangular_matrix.shape[1]]),2)
            for index in columns_to_delete:
                # If 0 on diagonal, insert 2 into solved_unknowns to represent erausre (cannot be string i.e. '?')
                solved_unknowns = np.insert(solved_unknowns, index, 2)
        except:
            # If not enough equations, use gauss jordan elimination to get reduced row echelon form.
            # Now, we can remove columns of resulting matrix which are not only 1 on the diagonal
            print("Didn't run")
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
                if solved_unknowns[0] == 2:
                    decoded_codeword += '?'
                else:
                    decoded_codeword += [solved_unknowns[0]]
                solved_unknowns = np.delete(solved_unknowns, 0)
        return decoded_codeword

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

    def new_message_pass_decode(self, binary_sequence, max_its):
        # Matrices that contain all message updates - not very memory efficient for larger codes
        Mcv = np.zeros_like(self.parity_check)
        Mvc = np.zeros(self.n)

        # Initialise v->c messages as the channel output for each variable v
        Mvc = np.array(binary_sequence)
        errors = []

        for it in range(max_its):
            if 0 in Mvc:
                # Update Mcv first
                for index,row in enumerate(self.parity_check):
                    variable_nodes = np.array(np.nonzero(row)[0])
                    for variable_node in variable_nodes:
                        other_variable_nodes = list(variable_nodes.copy())
                        other_variable_nodes.remove(variable_node)
                        Mcv[index][variable_node] = self.compute_cv_message(Mvc[other_variable_nodes])

                # Then update Mvc
                for index,column in enumerate(self.parity_check.T):
                    if Mvc[index] == 0:
                        check_nodes = np.array(np.nonzero(column)[0])
                        for check_node in check_nodes:
                            other_check_nodes = list(check_nodes.copy())
                            other_check_nodes.remove(check_node)
                            Mvc[index] = self.compute_vc_message(Mcv.T[index][other_check_nodes])
            errors.append(self.n - np.count_nonzero(Mvc))
        # Convert decoded values back to binary list, leaving '?' if erasures still present
        decoded_values = list(map(lambda x:'?' if x==0 else x, Mvc))
        return list(map(lambda x: 0 if x==-1 else x, decoded_values)), errors


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
        Mvc = np.zeros(self.n)

        # Initialise v->c messages as the channel output for each variable v
        for variable_node,check_node_list in variable_lookup.items():
            for check in check_node_list:
                Mvc[variable_node] = binary_sequence[variable_node]
        # print(Mvc)
        errors = []

        # Iterate for a fixed number of iterations, specified by max_its, skip message updates if no erasures in messages
        for it in range(max_its):
            if 0 in Mvc:
                # Update Mcv first:
                for check_node, variable_node_list in check_lookup.items():
                    for variable_node in variable_node_list:
                        # Remove current variable node from inputs to compute_cv_message function
                        other_variable_nodes = variable_node_list.copy()
                        other_variable_nodes.remove(variable_node)
                        Mcv[check_node][variable_node] = self.compute_cv_message(Mvc[other_variable_nodes])

                # Then update Mvc
                for variable_node, check_node_list in variable_lookup.items():
                    for check_node in check_node_list:
                        # If current Mvc message is already non-erasure, do nothing
                        # Otherwise compute v->c message
                        if Mvc[variable_node] == 0:
                            other_check_nodes = check_node_list.copy()
                            other_check_nodes.remove(check_node)
                            Mvc[variable_node] = self.compute_vc_message(Mcv.T[variable_node][other_check_nodes])
            errors.append(self.n - np.count_nonzero(Mvc))
        # Convert decoded values back to binary list, leaving '?' if erasures still present
        decoded_values = list(map(lambda x:'?' if x==0 else x, Mvc))
        return list(map(lambda x: 0 if x==-1 else x, decoded_values)), errors

def run_simulation(parameter_set):

    sim_BEC = BEC(parameter_set['BEC'])

    num_tests = parameter_set['num_tests']
    iterations = parameter_set['iterations']
    n = parameter_set['n']
    dv = parameter_set['dv']
    dc = parameter_set['dc']
    optimal = parameter_set['optimal']

    message_passing_block_errors = 0
    message_passing_bit_errors= 0
    if optimal:
        optimal_decoding_block_errors = 0
        optimal_decoding_bit_errors = 0

    average_errors = []
    skip = False

    for i in range(num_tests):
        parity_check = generate_random_parity_check(n,dv,dc)

        if len(parity_check) == 1:
            skip = True
            break
        else:
            skip = False

        LDPC = regular_LDPC_code(parity_check)
        codeword = np.zeros(LDPC.n)

        # print('Generated Codeword:', LDPC.encode(test_sequence))
        channel_output = sim_BEC.transmit(codeword)

        decoded_codeword, errors = LDPC.message_pass_decode(channel_output, iterations+1)
        if len(average_errors) == 0:
            average_errors += errors
            average_errors = np.array(average_errors)
        else:
            average_errors = average_errors*i + errors
            average_errors = average_errors / (i+1)

        if '?' in decoded_codeword:
            message_passing_block_errors += 1
        if optimal:
            optimal_decoded_codeword = LDPC.optimal_decode(channel_output)
            if '?' in optimal_decoded_codeword or -1 in optimal_decoded_codeword:
                optimal_decoding_block_errors += 1

            optimal_decoding_bit_errors += optimal_decoded_codeword.count('?')
        message_passing_bit_errors += decoded_codeword.count('?')
    if not skip:
        average_errors = average_errors / LDPC.n

        filename = 'regular_code'
        filename += '_BEC=' + str(sim_BEC.erasure_prob)
        filename += '_n=' + str(LDPC.n)
        filename += '_k=' + str(LDPC.k)
        filename += '_dv=' + str(LDPC.dv)
        filename += '_dc=' + str(LDPC.dc)
        filename += '_it=' + str(iterations)
        filename += '_num=' + str(num_tests)
        filename += '.csv'

        with open('./simulation_data/'+filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for error_at_iteration in average_errors:
                writer.writerow([error_at_iteration])
            writer.writerow(['Message passing block-wise error percentage', 100*message_passing_block_errors/num_tests])
            writer.writerow(['Message passing bit-wise error percentage', 100*message_passing_bit_errors/(num_tests*LDPC.n)])
            if optimal:
                writer.writerow(['Optimal decoding block-wise error percentage', 100*optimal_decoding_block_errors/num_tests])
                writer.writerow(['Optimal decoding bit-wise error percentage', 100*optimal_decoding_bit_errors/(num_tests*LDPC.n)])

        print('Message passing block-wise error percentage', 100*message_passing_block_errors/num_tests)
        if optimal:
            print('Optimal decoding block-wise error percentage', 100*optimal_decoding_block_errors/num_tests)

        print('Message passing bit-wise error percentage', 100*message_passing_bit_errors/(num_tests*LDPC.n))
        if optimal:
            print('Optimal decoding bit-wise error percentage', 100*optimal_decoding_bit_errors/(num_tests*LDPC.n)) 
        print(multiprocessing.current_process().name, ' done')

def run_simulation_fixed_ldpc(parameter_set):
    sim_BEC = BEC(parameter_set['BEC'])

    num_tests = parameter_set['num_tests']
    iterations = parameter_set['iterations']
    n = parameter_set['n']
    dv = parameter_set['dv']
    dc = parameter_set['dc']
    optimal = parameter_set['optimal']
    filenumber = parameter_set['filenumber']

    message_passing_block_errors = 0
    message_passing_bit_errors= 0
    if optimal:
        optimal_decoding_block_errors = 0
        optimal_decoding_bit_errors = 0

    average_errors = []
    parity_check = [0]
    while len(parity_check)==1:
        parity_check = generate_random_parity_check(n,dv,dc)

    LDPC = regular_LDPC_code(parity_check)
    codeword = np.zeros(LDPC.n)

    for i in range(num_tests):

        channel_output = sim_BEC.transmit(codeword)

        decoded_codeword, errors = LDPC.message_pass_decode(channel_output, iterations+1)

        if len(average_errors) == 0:
            average_errors += errors
            average_errors = np.array(average_errors)
        else:
            average_errors = average_errors*i + errors
            average_errors = average_errors / (i+1)

        if '?' in decoded_codeword:
            message_passing_block_errors += 1
        if optimal:
            optimal_decoded_codeword = LDPC.optimal_decode(channel_output)
            if '?' in optimal_decoded_codeword or -1 in optimal_decoded_codeword:
                optimal_decoding_block_errors += 1

            optimal_decoding_bit_errors += optimal_decoded_codeword.count('?')
        message_passing_bit_errors += decoded_codeword.count('?')

    average_errors = average_errors / LDPC.n

    filename = 'regular_code'
    filename += '_code_number=' + str(filenumber)
    filename += '_BEC=' + str(sim_BEC.erasure_prob)
    filename += '_n=' + str(LDPC.n)
    filename += '_k=' + str(LDPC.k)
    filename += '_dv=' + str(LDPC.dv)
    filename += '_dc=' + str(LDPC.dc)
    filename += '_it=' + str(iterations)
    filename += '_num=' + str(num_tests)
    filename += '.csv'

    with open('./simulation_data/'+filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for error_at_iteration in average_errors:
            writer.writerow([error_at_iteration])
        writer.writerow(['Message passing block-wise error percentage', 100*message_passing_block_errors/num_tests])
        writer.writerow(['Message passing bit-wise error percentage', 100*message_passing_bit_errors/(num_tests*LDPC.n)])
        if optimal:
            writer.writerow(['Optimal decoding block-wise error percentage', 100*optimal_decoding_block_errors/num_tests])
            writer.writerow(['Optimal decoding bit-wise error percentage', 100*optimal_decoding_bit_errors/(num_tests*LDPC.n)])

    print('Message passing block-wise error percentage', 100*message_passing_block_errors/num_tests)
    if optimal:
        print('Optimal decoding block-wise error percentage', 100*optimal_decoding_block_errors/num_tests)

    print('Message passing bit-wise error percentage', 100*message_passing_bit_errors/(num_tests*LDPC.n))
    if optimal:
        print('Optimal decoding bit-wise error percentage', 100*optimal_decoding_bit_errors/(num_tests*LDPC.n)) 
    print(multiprocessing.current_process().name, ' done')

ns = [50, 100, 200, 500, 1000]
# parameters = [{'BEC': 0.5, 'num_tests':10000, 'iterations':50, 'n':50, 'dv':3, 'dc':6}]
parameters = []
# for n in ns:
#     dictionary = {'BEC': 0.42, 'num_tests':10000, 'iterations':50, 'n':n, 'dv':3, 'dc':6, 'optimal':True}
#     parameters.append(dictionary)
#     dictionary = {'BEC': 0.43, 'num_tests':10000, 'iterations':50, 'n':n, 'dv':3, 'dc':6, 'optimal':True}
#     parameters.append(dictionary)

for i in range(32):
    dictionary = {'BEC': 0.42, 'num_tests':10000, 'iterations':50, 'n':50, 'dv':3, 'dc':6, 'optimal':False, 'filenumber':i}
    parameters.append(dictionary)

if __name__ == '__main__':
    pool = multiprocessing.Pool()
    for i in range(len(parameters)):
        pool.apply_async(run_simulation_fixed_ldpc, args=(parameters[i],))
    pool.close()
    pool.join()

# run_simulation_fixed_ldpc({'BEC': 0.42, 'num_tests':1000, 'iterations':50, 'n':100, 'dv':3, 'dc':6, 'optimal':False, 'filenumber':i})
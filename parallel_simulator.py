import numpy as np
from numpy.core.fromnumeric import var
from pyldpc.code import coding_matrix_systematic, coding_matrix
from pyldpc.utils import gausselimination, gaussjordan
from channels import BEC
from random_code_generator_python import generate_random_parity_check, generate_random_parity_check_no_checks
import matplotlib.pyplot as plt
import csv
import multiprocessing
import ctypes as ct
import os
import sys
from datetime import datetime

base_directory = '/home/rory/Documents/iib_project_ldpc_codes/'

def write_file(filename, errors, message_passing_block_error, message_passing_bit_error, optimal_block_error=None, optimal_bit_error=None):
    with open(base_directory + 'simulation_data/'+filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for error_at_iteration in errors:
            writer.writerow([error_at_iteration])
        writer.writerow(['Message passing block-wise error', message_passing_block_error])
        writer.writerow(['Message passing bit-wise error', message_passing_bit_error])
        if optimal_block_error:
            writer.writerow(['Optimal decoding block-wise error', optimal_block_error])
        if optimal_bit_error:
            writer.writerow(['Optimal decoding bit-wise error', optimal_bit_error])

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
        be found, returns input sequence.'''
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
            decoded_values = list(map(lambda x:'?' if x==0 else x, binary_sequence))
            return list(map(lambda x: 0 if x==-1 else x, decoded_values))
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
            decoded_values = list(map(lambda x:'?' if x==0 else x, binary_sequence))
            return list(map(lambda x: 0 if x==-1 else x, decoded_values))

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
                    decoded_codeword += [int(solved_unknowns[0])]
                solved_unknowns = np.delete(solved_unknowns, 0)
        return decoded_codeword

    def message_pass_decode(self, binary_sequence, max_its, check_lookup, variable_lookup):
        # Generate lookup tables for connected variable and check nodes ie dictionary 
        # with each check as key and connected variables list as value

        # check_lookup = []
        # for row in self.parity_check:
        #     check_list = list(np.nonzero(row==1)[0])
        #     check_lookup.append(check_list)
        # print(variable_lookup)
        # variable_lookup = []
        # for col in self.parity_check.T:
        #     var_list = list(np.nonzero(col==1)[0])
        #     variable_lookup.append(var_list)

        # print(variable_lookup)
        
        # Prepare variables for passing to C library
        check_lookup = np.array(check_lookup, dtype='int32').flatten()
        variable_lookup = np.array(variable_lookup, dtype='int32').flatten()
        errors = np.zeros(max_its, dtype='int32')
        initial_error_count = len(np.nonzero(binary_sequence==2)[0])
        binary_sequence = np.array(binary_sequence, dtype='int32')


        binary_sequence_p = binary_sequence.ctypes.data_as(ct.POINTER(ct.c_int))
        check_lookup_p = check_lookup.ctypes.data_as(ct.POINTER(ct.c_int))
        variable_lookup_p = variable_lookup.ctypes.data_as(ct.POINTER(ct.c_int))
        errors_p = errors.ctypes.data_as(ct.POINTER(ct.c_int))
        max_its_p = ct.c_int(max_its)
        n_p = ct.c_int(self.n)
        k_p = ct.c_int(self.k)
        dv_p = ct.c_int(self.dv)
        dc_p = ct.c_int(self.dc)

        c_message_pass = ct.CDLL(base_directory + 'message_passing.so')

        it = c_message_pass.message_passing(binary_sequence_p, max_its_p, variable_lookup_p, check_lookup_p, errors_p, n_p, k_p, dv_p, dc_p)
        errors = np.insert(errors, 0, initial_error_count)
        return binary_sequence, errors
        
        # decoded_values = list(map(lambda x:'?' if x==0 else x, Mvc))
        # return list(map(lambda x: 0 if x==-1 else x, decoded_values)), errors

def run_simulation(parameter_set):

    sim_BEC = BEC(parameter_set['BEC'])

    num_tests = parameter_set['num_tests']
    iterations = parameter_set['iterations']
    n = parameter_set['n']
    dv = parameter_set['dv']
    dc = parameter_set['dc']
    optimal = parameter_set['optimal']
    k = int(n*(dc-dv)/dc)

    message_passing_block_errors = 0
    message_passing_bit_errors= 0
    if optimal:
        optimal_decoding_block_errors = 0
        optimal_decoding_bit_errors = 0

    # average_errors = []
    error_counts = np.zeros(iterations+1)
    i = 0
    it = 0

    c_random_code = ct.CDLL(base_directory + 'random_code_generator.so')

    while i<num_tests and message_passing_block_errors<200:


        check_lookup = np.zeros(n*dv, dtype='int32')
        variable_lookup = np.zeros(n*dv, dtype='int32')
        parity_check = np.zeros(n*(n-k), dtype='bool')
        check_lookup_p = check_lookup.ctypes.data_as(ct.POINTER(ct.c_int))
        variable_lookup_p = variable_lookup.ctypes.data_as(ct.POINTER(ct.c_int))
        parity_check_p = parity_check.ctypes.data_as(ct.POINTER(ct.c_bool))

        it_p = ct.c_int(it)
        n_p = ct.c_int(n)
        dv_p = ct.c_int(dv)
        dc_p = ct.c_int(dc)

        success = c_random_code.generate_random_code(n_p, dv_p, dc_p, variable_lookup_p, check_lookup_p, parity_check_p, it_p)
        # print('Success', success)
        parity_check = np.reshape(parity_check, (n-k,n))
        # print(parity_check)

        # parity_check = generate_random_parity_check_no_checks(n,dv,dc)
        # if len(parity_check) == 1:
        #     continue
        print(i)

        LDPC = regular_LDPC_code(parity_check)
        codeword = np.zeros(LDPC.n)

        channel_output = sim_BEC.new_transmit(codeword)
        decoded_codeword, errors = LDPC.message_pass_decode(channel_output, iterations, check_lookup, variable_lookup)
        # decoded_codeword, errors = message_pass_decode(channel_output, LDPC.parity_check, iterations)

        error_counts += errors

        if errors[-1] != 0:
            message_passing_block_errors += 1
        if optimal:
            optimal_decoded_codeword = LDPC.optimal_decode(channel_output)
            if '?' in optimal_decoded_codeword or -1 in optimal_decoded_codeword:
                optimal_decoding_block_errors += 1

            optimal_decoding_bit_errors += optimal_decoded_codeword.count('?')
        message_passing_bit_errors += errors[-1]

        i += 1

    num_tests = i

    average_errors = error_counts / (LDPC.n*num_tests)

    filename = 'regular_code'
    filename += '_BEC=' + str(sim_BEC.erasure_prob)
    filename += '_n=' + str(LDPC.n)
    filename += '_k=' + str(LDPC.k)
    filename += '_dv=' + str(LDPC.dv)
    filename += '_dc=' + str(LDPC.dc)
    filename += '_it=' + str(iterations)
    filename += '_num=' + str(num_tests)
    filename += '_time=' + datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
    filename += '.csv'

    print("Printing to file")

    if optimal:
        write_file(filename, average_errors, message_passing_block_errors/num_tests, message_passing_bit_errors/(num_tests*LDPC.n), optimal_decoding_block_errors/num_tests, optimal_decoding_bit_errors/(num_tests*LDPC.n))
    else:
        write_file(filename, average_errors, message_passing_block_errors/num_tests, message_passing_bit_errors/(num_tests*LDPC.n))

    print(multiprocessing.current_process().name, ' done')

def run_simulation_fixed_ldpc(parameter_set):
    '''Run a simulation with a fixed code (for concentration plots)'''
    # Parse values from parameter_set dictionary
    BEC_value = parameter_set['BEC']
    num_tests = parameter_set['num_tests']
    iterations = parameter_set['iterations']
    n = parameter_set['n']
    dv = parameter_set['dv']
    dc = parameter_set['dc']
    optimal = parameter_set['optimal']
    filenumber = parameter_set['filenumber']

    # Check for existing numpy binary file containing parity check matrix
    filename = 'code_no_'+str(filenumber)+'_n_'+str(n)+'_dv_'+str(dv)+'_dc_'+str(dc)+'.npy'
    existing_code_file = False  
    for check_filename in os.listdir(base_directory+'parity_checks/'):
        if check_filename == filename:
            existing_code_file = True
            print('Found existing code')
            with open(base_directory+'parity_checks/' + filename, 'rb') as f:
                parity_check = np.load(f)

    # If no file found, generate new parity check and save it
    if not existing_code_file:
        parity_check = [0]
        while len(parity_check)==1:
            parity_check = generate_random_parity_check_no_checks(n,dv,dc)

        with open(base_directory+'parity_checks/' + filename, 'wb') as f:
            np.save(f, parity_check)

    parity_check = np.array(parity_check, dtype='int8')
    LDPC = regular_LDPC_code(parity_check)
    codeword = np.zeros(LDPC.n)

    message_passing_block_errors = 0
    message_passing_bit_errors= 0
    if optimal:
        optimal_decoding_block_errors = 0
        optimal_decoding_bit_errors = 0
    sim_BEC = BEC(BEC_value)
    average_errors = np.array([])

    i = 0

    while i<num_tests and message_passing_block_errors<200:
        if i%1000==0:
            print(i)
        channel_output = sim_BEC.new_transmit(codeword)

        decoded_codeword, errors = LDPC.message_pass_decode(channel_output, iterations)

        if len(average_errors) == 0:
            average_errors = np.array(errors)
        else:
            average_errors = average_errors*i + errors
            average_errors = average_errors / (i+1)

        if errors[-1] != 0:
            message_passing_block_errors += 1

        if optimal:
            optimal_decoded_codeword = LDPC.optimal_decode(channel_output)
            if '?' in optimal_decoded_codeword or -1 in optimal_decoded_codeword:
                optimal_decoding_block_errors += 1

            optimal_decoding_bit_errors += optimal_decoded_codeword.count('?')
        message_passing_bit_errors += errors[-1]
        i += 1
    average_errors = average_errors / LDPC.n

    filename = 'regular_code'
    filename += '_code_number=' + str(filenumber)
    filename += '_BEC=' + str(sim_BEC.erasure_prob)
    filename += '_n=' + str(LDPC.n)
    filename += '_k=' + str(LDPC.k)
    filename += '_dv=' + str(LDPC.dv)
    filename += '_dc=' + str(LDPC.dc)
    filename += '_it=' + str(iterations)
    filename += '_num=' + str(i)
    filename += '_time=' + datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
    filename += '.csv'

    if optimal:
        write_file(filename, average_errors, message_passing_block_errors/i, message_passing_bit_errors/(i*LDPC.n), optimal_decoding_block_errors/i, optimal_decoding_bit_errors/(i*LDPC.n))
    else:
        write_file(filename, average_errors, message_passing_block_errors/i, message_passing_bit_errors/(i*LDPC.n))

    print(multiprocessing.current_process().name, ' done')

parameters = []

erasure_prob = float(sys.argv[1])
num_tests = int(sys.argv[2])
iterations = int(sys.argv[3])
n = int(sys.argv[4])
dv = int(sys.argv[5])
dc = int(sys.argv[6])

# dictionary = {'BEC': erasure_prob, 'num_tests':num_tests, 'iterations':iterations, 'n':n, 'dv':dv, 'dc':dc, 'optimal':False, 'filenumber':filenumber}
# run_simulation_fixed_ldpc(dictionary)

dictionary = {'BEC': erasure_prob, 'num_tests':num_tests, 'iterations':iterations, 'n':n, 'dv':dv, 'dc':dc, 'optimal':False}
run_simulation(dictionary)
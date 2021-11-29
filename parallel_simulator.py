import numpy as np
from numpy.core.fromnumeric import var
from pyldpc.code import coding_matrix_systematic, coding_matrix
from pyldpc.utils import gausselimination, gaussjordan
from scipy.sparse import base
from channels import BEC
from random_code_generator_python import generate_random_parity_check, generate_random_parity_check_no_checks
import matplotlib.pyplot as plt
import csv
import multiprocessing
import ctypes as ct
import os
import sys
from datetime import datetime
import galois
import scipy.io

base_directory = '/home/rory/Documents/iib_project_ldpc_codes/'

def write_optimal_file(filename, optimal_block_error, optimal_bit_error):
    with open(base_directory + 'simulation_data/'+filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Optimal decoding block-wise error', optimal_block_error])
        writer.writerow(['Optimal decoding bit-wise error', optimal_bit_error])

def write_message_passing_file(filename, errors, message_passing_block_error, message_passing_bit_error):
    with open(base_directory + 'simulation_data/'+filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for error_at_iteration in errors:
            writer.writerow([error_at_iteration])
        writer.writerow(['Message passing block-wise error', message_passing_block_error])
        writer.writerow(['Message passing bit-wise error', message_passing_bit_error])

def write_combined_file(filename, errors, message_passing_block_error, message_passing_bit_error, optimal_block_error, optimal_bit_error):
    with open(base_directory + 'simulation_data/'+filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for error_at_iteration in errors:
            writer.writerow([error_at_iteration])
        writer.writerow(['Message passing block-wise error', message_passing_block_error])
        writer.writerow(['Message passing bit-wise error', message_passing_bit_error])
        writer.writerow(['Optimal decoding block-wise error', optimal_block_error])
        writer.writerow(['Optimal decoding bit-wise error', optimal_bit_error])

class regular_LDPC_code():

    def __init__(self,parity_check, n, k, dv, dc):
        # self.generator_matrix = coding_matrix(np.array(parity_check, dtype='int'))
        self.parity_check = parity_check
        self.n = n
        self.k =  k
        self.dv = dv
        self.dc = dc
        self.rate = self.k / self.n

    def encode(self,binary_sequence):
        '''Use generator matrix to encode the input binary sequence'''
        assert len(binary_sequence) == self.k, "The input sequence must have the correct input length (equal to k)"
        return np.mod(np.dot(self.generator_matrix, binary_sequence),2)

    def optimal_decode(self, binary_sequence):
        '''Use Gaussian Elimination and back subsitution to solve for the erased bits,
        given the values of parity check equations for the non-erased bits. If solution cannot
        be found, returns input sequence.'''

        # First check for trivial case of no erasures or more erasures than we have rows in the parity check
        no_erasures = np.count_nonzero(binary_sequence==2)
        if no_erasures==0 or no_erasures>(self.n-self.k):
            print("Either no erasures, or too many to be able to solve")
            return binary_sequence

        binary_sequence = np.array(binary_sequence, dtype='int32')
        target = np.zeros(self.n-self.k, dtype='bool')
        remaining_parity_checks = np.zeros(no_erasures*(self.n-self.k), dtype='bool')
        binary_sequence_p = binary_sequence.ctypes.data_as(ct.POINTER(ct.c_int))
        parity_check_p = self.parity_check.ctypes.data_as(ct.POINTER(ct.c_bool))
        remaining_parity_checks_p = remaining_parity_checks.ctypes.data_as(ct.POINTER(ct.c_bool))
        target_p = target.ctypes.data_as(ct.POINTER(ct.c_bool))
        n_p = ct.c_int(self.n)
        dv_p = ct.c_int(self.dv)
        dc_p = ct.c_int(self.dc)


        c_ml_decoder = ct.CDLL(base_directory + 'ml_decoder.so')
        c_ml_decoder.ml_decode(binary_sequence_p, target_p, parity_check_p, remaining_parity_checks_p, n_p, dv_p, dc_p)
        target = np.array(target, dtype='int')

        remaining_parity_checks = np.array(remaining_parity_checks, dtype='int')
        remaining_parity_checks = np.reshape(remaining_parity_checks, (no_erasures, self.n-self.k)).T
        GF = galois.GF(2)
        remaining_parity_checks_reduced = GF(np.c_[remaining_parity_checks, target])
        remaining_parity_checks_reduced = remaining_parity_checks_reduced.row_reduce(no_erasures)
        unsolvable_unknown_indexes = []
        erasure_positions = np.nonzero(binary_sequence==2)[0]
        i = 0
        # Check if reduced partity check matrix has variables we cannot solve
        while np.count_nonzero(remaining_parity_checks_reduced[:no_erasures,:-1].diagonal()==1) != (no_erasures-len(unsolvable_unknown_indexes)) and i<1000:
            # Get index in reduced parity check matrix of first unsolvable bit
            first_unknown_index = np.nonzero(remaining_parity_checks_reduced[:,:-1].diagonal()!=1)[0][0]
            # Convert this into the index in the overall codeword
            first_unknown_bit = erasure_positions[first_unknown_index]
            erasure_positions = np.delete(erasure_positions, first_unknown_index)
            # Add to list of unsolvables
            unsolvable_unknown_indexes.append(first_unknown_bit)
            # Get indexes of check nodes to remove
            checks_to_remove = np.nonzero(remaining_parity_checks[:,first_unknown_index])[0]
            remaining_parity_checks = np.delete(remaining_parity_checks, checks_to_remove, axis=0)
            remaining_parity_checks = np.delete(remaining_parity_checks, first_unknown_index, axis=1)
            target = np.delete(target, checks_to_remove)
            remaining_parity_checks_reduced = GF(np.c_[remaining_parity_checks, target]).row_reduce(no_erasures-len(unsolvable_unknown_indexes))
            i += 1
            # print("Unsolvables: ", unsolvable_unknown_indexes)

        solved_unknowns = np.array(remaining_parity_checks_reduced[:no_erasures-len(unsolvable_unknown_indexes),-1])
        
        decoded_codeword = []
        for index,bit in enumerate(binary_sequence):
            if bit == 2:
                if index not in unsolvable_unknown_indexes:
                    decoded_codeword.append(solved_unknowns[0])
                    solved_unknowns = np.delete(solved_unknowns, 0)
                else:
                    # print(index)
                    decoded_codeword.append(2)
            else:
                decoded_codeword.append(bit)
            if bit == 0 and decoded_codeword[index]==2:
                print(index, bit, decoded_codeword[index])

        return np.array(decoded_codeword)

    def message_pass_decode(self, binary_sequence, max_its, check_lookup = None, variable_lookup = None):
        if check_lookup is None:
            check_lookup = []
            for row in self.parity_check:
                check_list = list(np.nonzero(row==1)[0])
                check_lookup.append(check_list)

        if variable_lookup is None:
            variable_lookup = []
            for col in self.parity_check.T:
                var_list = list(np.nonzero(col==1)[0])
                variable_lookup.append(var_list)
       
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

def run_simulation(parameter_set):

    sim_BEC = BEC(parameter_set['BEC'])
    num_tests = parameter_set['num_tests']
    iterations = parameter_set['iterations']
    n = parameter_set['n']
    dv = parameter_set['dv']
    dc = parameter_set['dc']
    optimal = parameter_set['optimal']
    message_passing = parameter_set['message_passing']
    seed = parameter_set['seed']
    k = int(n*(dc-dv)/dc)
    start_time = datetime.now()

    if message_passing:
        message_passing_block_errors = 0
        message_passing_bit_errors= 0
        block_error = message_passing_block_errors
    if optimal:
        optimal_decoding_block_errors = 0
        optimal_decoding_bit_errors = 0
        if not message_passing:
            block_error = optimal_decoding_block_errors

    message_passing_error_counts = np.zeros(iterations+1)
    i = 0
    it = 0

    c_random_code = ct.CDLL(base_directory + 'random_code_generator.so')

    while block_error<200 and i<num_tests and (datetime.now()-start_time).total_seconds() < 43000:

        # Prepare variables for C random code generator
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
        first_run = i==0
        first_run_p = ct.c_bool(first_run)
        seed_p = ct.c_int(seed)

        success = c_random_code.generate_random_code(n_p, dv_p, dc_p, variable_lookup_p, check_lookup_p, parity_check_p, it_p, first_run_p, seed_p)
        parity_check = np.reshape(parity_check, (n-k,n))
        # if i%1000==0:
        print(i)

        # Create LDPC code and transmit over channel
        LDPC = regular_LDPC_code(parity_check, n, k, dv, dc)
        codeword = np.zeros(LDPC.n)
        channel_output = sim_BEC.new_transmit(codeword)

        if message_passing:
            decoded_codeword, errors = LDPC.message_pass_decode(channel_output, iterations, check_lookup, variable_lookup)
            message_passing_error_counts += errors
            if errors[-1] != 0:
                message_passing_block_errors += 1
            message_passing_bit_errors += errors[-1]
            block_error = message_passing_block_errors

        if optimal:
            optimal_decoded_codeword = LDPC.optimal_decode(channel_output)
            optimal_error_count = np.count_nonzero(optimal_decoded_codeword==2)

            if optimal_error_count > 0:
                optimal_decoding_block_errors += 1
            optimal_decoding_bit_errors += optimal_error_count

            if not message_passing:
                block_error = optimal_decoding_block_errors

        i += 1

    num_tests = i

    average_errors = message_passing_error_counts / (LDPC.n*num_tests)

    filename = 'regular_code'
    filename += '_BEC=' + str(sim_BEC.erasure_prob)
    filename += '_n=' + str(LDPC.n)
    filename += '_k=' + str(LDPC.k)
    filename += '_dv=' + str(LDPC.dv)
    filename += '_dc=' + str(LDPC.dc)
    if message_passing:
        filename += '_it=' + str(iterations)
    filename += '_num=' + str(num_tests)
    filename += '_time=' + datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
    filename += '.csv'

    print("Printing to file")

    if optimal and message_passing:
        write_combined_file(filename, average_errors, message_passing_block_errors/num_tests, message_passing_bit_errors/(num_tests*LDPC.n), optimal_decoding_block_errors/num_tests, optimal_decoding_bit_errors/(num_tests*LDPC.n))
    elif optimal:
        write_optimal_file(filename, optimal_decoding_block_errors/num_tests, optimal_decoding_bit_errors/(num_tests*LDPC.n))
    elif message_passing:
        write_message_passing_file(filename, average_errors, message_passing_block_errors/num_tests, message_passing_bit_errors/(num_tests*LDPC.n))


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
    message_passing = parameter_set['message_passing']
    filenumber = parameter_set['filenumber']
    k = int(n*(dc-dv)/dc)
    start_time = datetime.now()

    # Check for existing numpy binary file containing parity check matrix
    code_filename = 'code_no_'+str(filenumber)+'_n_'+str(n)+'_dv_'+str(dv)+'_dc_'+str(dc)+'.npy'
    check_filename = 'check_code_no_'+str(filenumber)+'_n_'+str(n)+'_dv_'+str(dv)+'_dc_'+str(dc)+'.npy'
    variable_filename = 'variable_code_no_'+str(filenumber)+'_n_'+str(n)+'_dv_'+str(dv)+'_dc_'+str(dc)+'.npy'
    files_found = 0
    for filename in os.listdir(base_directory+'parity_checks/'):
        if filename == code_filename:
            print('Found existing code')
            files_found += 1
            with open(base_directory+'parity_checks/' + code_filename, 'rb') as f:
                parity_check = np.load(f)
        if filename == check_filename:
            files_found += 1
            print('Found existing check lookup')
            with open(base_directory+'parity_checks/' + check_filename, 'rb') as f:
                check_lookup = np.load(f)
        if filename == variable_filename:
            files_found += 1
            print('Found existing variable lookup')
            with open(base_directory+'parity_checks/' + variable_filename, 'rb') as f:
                variable_lookup = np.load(f)

    # If no file found, generate new parity check and save it
    if files_found != 3:
        it = 0
        c_random_code = ct.CDLL(base_directory + 'random_code_generator.so')

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
        parity_check = np.reshape(parity_check, (n-k,n))

        with open(base_directory+'parity_checks/' + code_filename, 'wb') as f:
            np.save(f, parity_check)
        with open(base_directory+'parity_checks/' + check_filename, 'wb') as f:
            np.save(f, check_lookup)
        with open(base_directory+'parity_checks/' + variable_filename, 'wb') as f:
            np.save(f, variable_lookup)

    LDPC = regular_LDPC_code(parity_check, n, k, dv, dc)
    codeword = np.zeros(LDPC.n)

    if message_passing:
        message_passing_block_errors = 0
        message_passing_bit_errors= 0
        block_error = message_passing_block_errors
    if optimal:
        optimal_decoding_block_errors = 0
        optimal_decoding_bit_errors = 0
        if not message_passing:
            block_error = optimal_decoding_block_errors

    message_passing_error_counts = np.zeros(iterations+1)
    sim_BEC = BEC(BEC_value)
    i = 0

    while block_error<200 and i<num_tests and (datetime.now()-start_time).total_seconds() < 42000:
        if i%100==0:
            print(i)
        channel_output = sim_BEC.new_transmit(codeword)

        if message_passing:
            decoded_codeword, errors = LDPC.message_pass_decode(channel_output, iterations, check_lookup, variable_lookup)

            message_passing_error_counts += errors
            if errors[-1] != 0:
                message_passing_block_errors += 1
            message_passing_bit_errors += errors[-1]
            block_error = message_passing_block_errors

        if optimal:
            optimal_decoded_codeword = LDPC.optimal_decode(channel_output)
            optimal_error_count = np.count_nonzero(optimal_decoded_codeword==2)

            if optimal_error_count > 0:
                optimal_decoding_block_errors += 1
            optimal_decoding_bit_errors += optimal_error_count

            if not message_passing:
                block_error = optimal_decoding_block_errors

        i += 1
    num_tests = i

    average_errors = message_passing_error_counts / (LDPC.n*num_tests)

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

    if optimal and message_passing:
        write_combined_file(filename, average_errors, message_passing_block_errors/num_tests, message_passing_bit_errors/(num_tests*LDPC.n), optimal_decoding_block_errors/num_tests, optimal_decoding_bit_errors/(num_tests*LDPC.n))
    elif optimal:
        write_optimal_file(filename, optimal_decoding_block_errors/num_tests, optimal_decoding_bit_errors/(num_tests*LDPC.n))
    elif message_passing:
        write_message_passing_file(filename, average_errors, message_passing_block_errors/num_tests, message_passing_bit_errors/(num_tests*LDPC.n))

parameters = []

erasure_prob = float(sys.argv[1])
num_tests = int(sys.argv[2])
iterations = int(sys.argv[3])
n = int(sys.argv[4])
dv = int(sys.argv[5])
dc = int(sys.argv[6])
mode = int(sys.argv[7])


if mode == 0:
    # Only message passing for random codes
    seed = int(sys.argv[8])
    dictionary = {'BEC': erasure_prob, 'num_tests':num_tests, 'iterations':iterations, 'n':n, 'dv':dv, 'dc':dc, 'optimal':False, 'message_passing': True, 'seed':seed}
    run_simulation(dictionary)
elif mode == 1:
    # Only ML decoder for random codes
    seed = int(sys.argv[8])
    dictionary = {'BEC': erasure_prob, 'num_tests':num_tests, 'iterations':iterations, 'n':n, 'dv':dv, 'dc':dc, 'optimal':True, 'message_passing': False, 'seed':seed}
    run_simulation(dictionary)
elif mode == 2:
    # Both ML decoder and message passing for random codes
    seed = int(sys.argv[8])
    dictionary = {'BEC': erasure_prob, 'num_tests':num_tests, 'iterations':iterations, 'n':n, 'dv':dv, 'dc':dc, 'optimal':True, 'message_passing': True, 'seed':seed}
    run_simulation(dictionary)
elif mode == 3:
    # Only message passing for fixed code
    filenumber = int(sys.argv[8])
    dictionary = {'BEC': erasure_prob, 'num_tests':num_tests, 'iterations':iterations, 'n':n, 'dv':dv, 'dc':dc, 'optimal':False, 'message_passing': True, 'filenumber':filenumber}
    run_simulation_fixed_ldpc(dictionary)
elif mode == 4:
    # Only ML decoder for fixed code
    filenumber = int(sys.argv[8])
    dictionary = {'BEC': erasure_prob, 'num_tests':num_tests, 'iterations':iterations, 'n':n, 'dv':dv, 'dc':dc, 'optimal':True,'message_passing': False, 'filenumber':filenumber}
    run_simulation_fixed_ldpc(dictionary)
elif mode == 5:
    # Both ML decoder and message passing for fixed code
    filenumber = int(sys.argv[8])
    dictionary = {'BEC': erasure_prob, 'num_tests':num_tests, 'iterations':iterations, 'n':n, 'dv':dv, 'dc':dc, 'optimal':True,'message_passing': True, 'filenumber':filenumber}
    run_simulation_fixed_ldpc(dictionary)
else:
    raise ValueError("Mode value must be in the range 0-5")


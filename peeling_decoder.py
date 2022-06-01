import numpy as np
from channels import BEC
import matplotlib.pyplot as plt
import ctypes as ct
import csv
from test_de_threshold import calc_threshold
import scipy.stats as stats
import math

base_directory = '/Users/rory/Documents/iib_project_ldpc_codes_new/iib_project_ldpc_codes/'

def calculate_threshold_y(threshold_erasure,dv,dc):
    previous_value = 0
    test_value = 1
    while np.abs(test_value - previous_value) > 1e-6:
        previous_value = test_value
        test_value = 1 - (1 - threshold_erasure * test_value**(dv-1))**(dc-1)
    return test_value
        
def generate_random_parity_check(n,k,dv,dc):
    c_random_code = ct.CDLL(base_directory + 'random_code_generator.so')
    sequence = np.arange(n*dv, dtype='int32')
    sequence_p = sequence.ctypes.data_as(ct.POINTER(ct.c_int))

    # Prepare variables for C random code generator
    check_lookup = np.zeros(n*dv, dtype='int32')
    variable_lookup = np.zeros(n*dv, dtype='int32')
    parity_check = np.zeros(n*(n-int(k)), dtype='bool')
    check_lookup_p = check_lookup.ctypes.data_as(ct.POINTER(ct.c_int))
    variable_lookup_p = variable_lookup.ctypes.data_as(ct.POINTER(ct.c_int))
    parity_check_p = parity_check.ctypes.data_as(ct.POINTER(ct.c_bool))
    it = 0
    seed = np.random.randint(1000)
    it_p = ct.c_int(it)
    n_p = ct.c_int(n)
    dv_p = ct.c_int(dv)
    dc_p = ct.c_int(dc)
    first_run = 0
    first_run_p = ct.c_bool(first_run)
    seed_p = ct.c_int(seed)

    success = c_random_code.generate_random_code(n_p, dv_p, dc_p, variable_lookup_p, check_lookup_p, sequence_p, parity_check_p, it_p, first_run_p, seed_p)
    
    parity_check = np.array(parity_check, dtype='int')
    return np.reshape(parity_check, (n-int(k),n))

def decode(received_message, n, dv, dc, parity_check):
    check_node_degrees = np.array([dc] * parity_check.shape[0])
    num_variable_nodes = n
    one_degree_evolution = []

    # Go through the received message, deleting items from parity check as you go
    for bit_index,bit in enumerate(received_message):
        # If bit is not erased, we need to find the corresponding check nodes and remove them
        if bit != 2:
            # Get the column of the non-erased bit
            for col_index,value in enumerate(parity_check[:,bit_index]):
                if value == 1:
                    check_node_degrees[col_index] -= 1
                    parity_check[col_index,bit_index] = 0
            num_variable_nodes -= 1

    while 1 in check_node_degrees:
        one_degree_indexes = [i for i in range(len(check_node_degrees)) if check_node_degrees[i] == 1]
        one_degree_evolution.append(len(one_degree_indexes))
        random_choice = np.random.choice(one_degree_indexes)
        for index in range(len(parity_check[random_choice,:])):
            if parity_check[random_choice,:][index] == 1:
                # We have found the variable node to be corrected
                # Now decrease the connected check node degrees
                column = parity_check[:,index]
                for check_node_index in range(len(column)):
                    if column[check_node_index] == 1:
                        check_node_degrees[check_node_index] -= 1
                parity_check[:,index] = 0
                break

    # Take care of last variable node 
    if not check_node_degrees.any():
        one_degree_evolution.append(0)

    return num_variable_nodes, one_degree_evolution

def calculate_alpha(threshold, dv, dc):
    y_crit = calculate_threshold_y(threshold, dv,dc)
    x_crit = threshold * (y_crit)**(dv-1)
    return threshold * np.sqrt( ((dv-1)/dv) * (1/x_crit - 1/y_crit)) 

n = 500

dv = 3
dc = 6
erasure_prob = 0.34
repeats = 100
results = []
max_length = 0
BEC = BEC(erasure_prob)
max_value = 0
min_value = 0

def f(erasure_prob,dv,dc,y):
    y = (1 - y/(erasure_prob*n))**(1/dv)
    return erasure_prob * y**(dv-1) * (y - 1 + (1-erasure_prob*y**(dv-1))**(dc-1))

threshold = calc_threshold(dv,dc)
ys = np.arange(int(n*erasure_prob))
drift = dv*n*f(erasure_prob,dv,dc,ys)[::-1]

y_crit = calculate_threshold_y(threshold, 3,6)
print(y_crit)
critical_point = int(n*threshold*(y_crit)**(dv))
critical_distribution = []
num_failures_at_critical = 0

sizes_at_failure = []

first_derivative_term = -n * dv * (dc-1) * threshold * y_crit**(2*dv-2) * (1-erasure_prob*(y_crit)**(dv-1))**(dc-2)
calculation = first_derivative_term * (erasure_prob - threshold)

second_derivative_term = (-2 * first_derivative_term / threshold) + n * dv * threshold * (dc-1) * (dc-2) * (y_crit)**(3*dv-3) * (1-erasure_prob*(y_crit)**(dv-1))**(dc-3)

second_derivative_term = n * dv * (dc-1) * (y_crit)**(2*dv-2) * (1-threshold*(y_crit)**(dv-1))**(dc-3) * (threshold * (dc-2) * y_crit**(dv-1) - 2*(1-threshold * y_crit**(dv-1)))
second_order_approx = calculation + second_derivative_term * 0.5*(erasure_prob - threshold)**2

# calculation = dv * n *(threshold-erasure_prob) * threshold * (y_crit)**(2*dv-2) * (1-erasure_prob*(y_crit)**(dv-1))**(dc-2)
print("Calculation", calculation)
print("Value", drift[int(critical_point)])

print("Second order",second_order_approx )

plt.plot(ys,drift)
plt.vlines(critical_point, 0, max(drift))
plt.show()


for i in range(repeats):
    # print(i)
    if i % 10 ==0:
        print(i)
    parity_check = generate_random_parity_check(n,n/2,dv,dc)
    received = BEC.new_transmit(np.zeros(n))
    num_variable_nodes, one_degree_evolution = decode(received, n, dv, dc, parity_check)
    sizes_at_failure.append(num_variable_nodes+1-len(one_degree_evolution))
    one_degree_evolution += [np.NaN]*(num_variable_nodes+1-len(one_degree_evolution))
    results.append(one_degree_evolution[::-1])
    if len(one_degree_evolution) > max_length:
        max_length = len(one_degree_evolution)

    one_degree_evolution = one_degree_evolution[::-1]

    plt.plot(one_degree_evolution, '--')
    result = []
    for j in range(len(ys)):
        if j == len(one_degree_evolution):
            break
        # if one_degree_evolution[j] == nan:
        #     result.append(np.nan)
        #     print("Appending nan")
        result.append(one_degree_evolution[j]-drift[j])
        if result[-1] > max_value:
            max_value = result[-1]
        if result[-1] < min_value:
            min_value = result[-1]
        if j == critical_point:
            if math.isnan(one_degree_evolution[j]-drift[j]) or one_degree_evolution[j]-drift[j] < 0:
                num_failures_at_critical += 1
            else:
                critical_distribution.append(one_degree_evolution[j])

    # plt.plot(result, '--')
    # plt.plot(one_degree_evolution, '--')
    # plt.plot(one_degree_evolution, '--', label=str(i))

# plt.plot(ys,drift)
ys = np.linspace(0,n*erasure_prob,1000)
plt.plot(ys, dv*n*f(erasure_prob,dv,dc,ys)[::-1])
# results.append(dv*n*f(erasure_prob,dv,dc,ys)[::-1])

# results = np.array(results)
# with open(base_directory + 'finite_length_scaling_report/peeling_decoder_path_simulation_n='+str(n)+'.csv', 'w') as csv_file:
#     csv_writer = csv.writer(csv_file)
    
#     for index in range(max_length):
#         data = [index]
#         for result in results[:-1]:
#             try:
#                 data.append(result[index])
#             except:
#                 data.append(np.nan)

#         csv_writer.writerow(data)

# with open(base_directory + 'finite_length_scaling_report/peeling_decoder_path_expected_n='+str(n)+'.csv', 'w') as csv_file:
#     csv_writer = csv.writer(csv_file)

#     for index,result in enumerate(results[-1]):
#         csv_writer.writerow([ys[index],result])

threshold = calc_threshold(3,6)
print(threshold)
y_crit = calculate_threshold_y(threshold, 3,6)
print(y_crit)
print(n*threshold*(y_crit)**(dv))
# plt.vlines(n*threshold*(y_crit)**(dv), min_value, max_value)



# print(calculate_threshold_y(threshold, 3,6))
# y = (1 - y*dv/(erasure_prob*n*dv))**(1/dv)

plt.legend()
plt.show()

print(len(critical_distribution))
# critical_distribution += [0] * num_failures_at_critical
print(len(critical_distribution))
print(num_failures_at_critical)


plt.hist(critical_distribution, bins=40)
alpha = calculate_alpha(threshold, dv, dc)
print(alpha)
x = np.linspace(min(critical_distribution),max(critical_distribution),1000)

variance = n * dv**2 * (alpha * threshold * (dc-1) * (y_crit)**(2*dv-2) * (1 - threshold * y_crit**(dv-1))**(dc-2))**2
print(variance)

plt.plot(x, (len(critical_distribution) + num_failures_at_critical) * stats.norm.pdf(x, loc=calculation, scale=np.sqrt(variance)))
plt.show()

plt.hist(sizes_at_failure)
plt.show()
import os
import numpy as np
import csv
data_directory = '/home/rory/Documents/iib_project_ldpc_codes/hpc_simulation_data_new/'

seen_files = []

def get_parameters(filename):
    # Get simulation parameters from filename
    simulation_parameters = filename.split('_')
    if 'combined.csv' in simulation_parameters:
        return {'BEC':'0.0', 'number':'1'}
    simulation_parameters_dict = {}
    for simulation_parameter in simulation_parameters:
        if '=' in simulation_parameter and '.csv' not in simulation_parameter:
            simulation_parameters_dict[simulation_parameter.split('=')[0]] = simulation_parameter.split('=')[1]
        elif '.csv' in simulation_parameter:
            simulation_parameters_dict[simulation_parameter.split('=')[0]] = simulation_parameter.split('=')[1].split('.')[0]
    return simulation_parameters_dict


    # print(get_parameters(filename))
BECs = [0.3, 0.32]

for code_no in range(1,11):
    for BEC in BECs:
        errors = np.zeros(200, dtype='int')
        file_count = 0
        num_count = 0
        for filename in os.listdir(data_directory):
            parameters = get_parameters(filename)
            if float(parameters['BEC']) == BEC and int(parameters['number']) == code_no:
                specific_error = []
                file_count += 1
                num_count += int(parameters['num'])
                with open(data_directory + filename) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    for row in csv_reader:
                        if len(row) != 1:
                            break
                        specific_error.append(round(float(row[0])*int(parameters['num'])*int(parameters['n'])))
                errors += np.array(specific_error)
        errors = errors / (num_count*int(parameters['n']))
        print('Code_no: ', code_no)
        print('BEC: ', BEC)
        print(errors)
        new_filename = 'regular_code_code_number=' + str(code_no) + '_BEC=' + str(BEC) + '_n='+parameters['n'] + '_k=' + parameters['k'] + '_dv='+parameters['dv'] + '_dc='+parameters['dc']+'_it='+parameters['it']+'_combined.csv'
        print(new_filename)
        
        with open(data_directory+new_filename, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            for error in errors:
                writer.writerow([error])
                print(error)

# with open('./simulation_data/'+filename, 'w', newline='') as csvfile:
#         writer = csv.writer(csvfile)
#         for error_at_iteration in errors:
#             writer.writerow([error_at_iteration])
#         writer.writerow(['Message passing block-wise error', message_passing_block_error])
#         writer.writerow(['Message passing bit-wise error', message_passing_bit_error])
#         if optimal_block_error:
#             writer.writerow(['Optimal decoding block-wise error', optimal_block_error])
#         if optimal_bit_error:
#             writer.writerow(['Optimal decoding bit-wise error', optimal_bit_error])
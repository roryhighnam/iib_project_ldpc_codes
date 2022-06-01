import os
import numpy as np
import csv
# data_directory = '/home/rory/Documents/iib_project_ldpc_codes/hpc_simulation_data_fixed/'
# return_combined_directory = '/home/rory/Documents/iib_project_ldpc_codes/hpc_simulation_data_fixed_mp_512/50_it/'

data_directory = '/Users/Rory/Documents/iib_project_ldpc_codes_new/iib_project_ldpc_codes/simulation_data_random_ensemble_report/'
return_combined_directory = '/Users/Rory/Documents/iib_project_ldpc_codes_new/iib_project_ldpc_codes/random_ensemble_report_combined/'

# if data_directory == return_combined_directory:
#     raise ValueError('Return directory must be different!')

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
    print(simulation_parameters_dict)
    return simulation_parameters_dict


    # print(get_parameters(filename))

def concentration_combine(BECs):
    for code_no in range(1,11):
        for BEC in BECs:
            errors = np.zeros(501, dtype='int')
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
            
            with open(return_combined_directory+new_filename, 'w', newline='') as csv_file:
                writer = csv.writer(csv_file)
                for error in errors:
                    writer.writerow([error])
                    print(error)

def ensemble_combine(BECs):
    for BEC in BECs:
        errors = np.zeros(202, dtype='int')
        file_count = 0
        num_count = 0
        for filename in os.listdir(data_directory):
            parameters = get_parameters(filename)
            print(parameters)
            if 'BEC' in parameters and 'id' in parameters and float(parameters['BEC']) == BEC:
                print(filename)
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
        # print('Code_no: ', code_no)
        print('BEC: ', BEC)
        print(errors)
        new_filename = 'regular_code' + '_BEC=' + str(BEC) + '_n='+parameters['n'] + '_k=' + parameters['k'] + '_dv='+parameters['dv'] + '_dc='+parameters['dc']+'_it='+parameters['it']+'_num='+str(num_count)+'.csv'
        print(new_filename)
        
        with open(return_combined_directory+new_filename, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            for error in errors:
                writer.writerow([error])
                print(error)

def optimal_concentration_combine(BECs):
    for code_no in range(1,11):
        for BEC in BECs:
            errors = 0
            file_count = 0
            num_count = 0
            for filename in os.listdir(data_directory):
                parameters = get_parameters(filename)
                if float(parameters['BEC']) == BEC and int(parameters['number']) == code_no:
                    specific_error = 0
                    file_count += 1
                    num_count += int(parameters['num'])
                    with open(data_directory + filename) as csv_file:
                        csv_reader = csv.reader(csv_file, delimiter=',')
                        for row in csv_reader:
                            if csv_reader.line_num != 1:
                                specific_error = round(float(row[1])*int(parameters['num'])*int(parameters['n']))
                    errors += specific_error
            errors = errors / (num_count*int(parameters['n']))
            print('Code_no: ', code_no)
            print('BEC: ', BEC)
            print(errors)
            new_filename = 'regular_code_code_number=' + str(code_no) + '_BEC=' + str(BEC) + '_n='+parameters['n'] + '_k=' + parameters['k'] + '_dv='+parameters['dv'] + '_dc='+parameters['dc']+'_it='+parameters['it']+'_num='+str(num_count)+'_combined.csv'
            print(new_filename)
            
            with open(return_combined_directory+new_filename, 'w', newline='') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow([errors])

def mp_concentration_combine(BECs, iteration):
    for code_no in range(1,11):
        for BEC in BECs:
            errors = 0
            file_count = 0
            num_count = 0
            for filename in os.listdir(data_directory):
                parameters = get_parameters(filename)
                if float(parameters['BEC']) == BEC and int(parameters['number']) == code_no:
                    specific_error = 0
                    file_count += 1
                    num_count += int(parameters['num'])
                    with open(data_directory + filename) as csv_file:
                        csv_reader = csv.reader(csv_file, delimiter=',')
                        for row in csv_reader:
                            if csv_reader.line_num == iteration:
                                specific_error = round(float(row[0])*int(parameters['num'])*int(parameters['n']))
                                break
                    errors += specific_error
            errors = errors / (num_count*int(parameters['n']))
            print('Code_no: ', code_no)
            print('BEC: ', BEC)
            print(errors)
            new_filename = 'regular_code_code_number=' + str(code_no) + '_BEC=' + str(BEC) + '_n='+parameters['n'] + '_k=' + parameters['k'] + '_dv='+parameters['dv'] + '_dc='+parameters['dc']+'_it='+str(iteration)+'_num='+str(num_count)+'_combined.csv'
            print(new_filename)
            
            with open(return_combined_directory+new_filename, 'w', newline='') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow([errors])                

def optimal_and_message_passing_concentration_combine(BECs):
    for code_no in range(11,31):
        print(code_no)
        for BEC in BECs:
            print(BEC)
            ml_errors = 0
            message_passing_errors = 0
            file_count = 0
            num_count = 0
            for filename in os.listdir(data_directory):
                parameters = get_parameters(filename)
                if float(parameters['BEC']) == BEC and int(parameters['number']) == code_no:
                    print(filename)
                    specific_ml_errors = 0
                    specific_message_passing_errors = 0
                    file_count += 1
                    num_count += int(parameters['num'])
                    with open(data_directory+filename) as csv_file:
                        csv_reader = csv.reader(csv_file, delimiter=',')
                        for row in csv_reader:
                            if csv_reader.line_num == 203:
                                specific_message_passing_errors = round(float(row[1])*int(parameters['num'])*int(parameters['n']))
                            if csv_reader.line_num == 205:
                                specific_ml_errors = round(float(row[1])*int(parameters['num'])*int(parameters['n']))

                    ml_errors += specific_ml_errors
                    message_passing_errors += specific_message_passing_errors

            ml_errors = ml_errors / (num_count*int(parameters['n']))
            message_passing_errors = message_passing_errors / (num_count*int(parameters['n']))
            new_filename = 'regular_code_code_number=' + str(code_no) + '_BEC=' + str(BEC) + '_n='+parameters['n'] + '_k=' + parameters['k'] + '_dv='+parameters['dv'] + '_dc='+parameters['dc']+'_it='+parameters['it']+'_num='+str(num_count)+'_combined.csv'
            with open(return_combined_directory+new_filename, 'w', newline='') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow([ml_errors])
                writer.writerow([message_passing_errors])
# def both_concentration_combine(BECs):

# ensemble_combine([0.42])
# optimal_concentration_combine([0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5])
# optimal_and_message_passing_concentration_combine([0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5])
ensemble_combine([0.4])
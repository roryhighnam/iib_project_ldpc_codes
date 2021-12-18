import numpy as np
import os

# base_directory = 'home/rjh248/iib_project_ldpc_codes/'
base_directory = '/home/rory/Documents/iib_project_ldpc_codes/'


# filenames = []
# for i in range(11):
#     filenames.append('code_no_'+str(i)+'_n_512_dv_3_dc_6.npy')

filenames = ['code_no_4_n_512_dv_3_dc_6.npy']

for filename in filenames:
    for test_filename in os.listdir(base_directory+'hpc_parity_checks/'):
        if filename == test_filename:
            print(filename)
            existing_code_file = True
            print('Found existing code')
            with open(base_directory+'hpc_parity_checks/' + filename, 'rb') as f:
                parity_check = np.load(f)
            print(parity_check)
            check_lookup = []
            for row in parity_check:
                check_list = list(np.nonzero(row==1)[0])
                check_lookup.append(check_list)

            variable_lookup = []
            for col in parity_check.T:
                var_list = list(np.nonzero(col==1)[0])
                variable_lookup.append(var_list)
            with open(base_directory+'hpc_parity_checks/check_' + filename, 'wb') as f:
                np.save(f, check_lookup)
            with open(base_directory+'hpc_parity_checks/variable_' + filename, 'wb') as f:
                np.save(f, variable_lookup)

    for i in range(len(variable_lookup)):
        for j in range(len(variable_lookup)):
            if i!= j:
                identical = True
                for l in range(3):
                    if variable_lookup[i][l] not in variable_lookup[j]:
                        identical = False
                if identical:
                    print(variable_lookup[i])
                    print(variable_lookup[j])
                    print(i,j)
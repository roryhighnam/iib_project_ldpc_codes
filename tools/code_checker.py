import numpy as np
from numpy.core.fromnumeric import var


code_directory = '/home/rory/Documents/iib_project_ldpc_codes/hpc_parity_checks_1000/'

for code in range(1,11):
    print(code)
    with open(code_directory+'code_no_'+str(code)+'_n_1000_dv_3_dc_6.npy', 'rb') as f:
        parity_check = np.load(f)
    with open(code_directory+'check_code_no_'+str(code)+'_n_1000_dv_3_dc_6.npy', 'rb') as f:
        stored_check_lookup = np.load(f)
    with open(code_directory+'variable_code_no_'+str(code)+'_n_1000_dv_3_dc_6.npy', 'rb') as f:
        stored_variable_lookup = np.load(f)

    check_lookup = []
    for row in parity_check:
        check_list = list(np.nonzero(row==1)[0])
        check_lookup.append(check_list)
    check_lookup = np.array(check_lookup).flatten()

    variable_lookup = []
    for col in parity_check.T:
        var_list = list(np.nonzero(col==1)[0])
        variable_lookup.append(var_list)
    variable_lookup = np.array(variable_lookup).flatten()

    zero_variable_checks = stored_variable_lookup[:3]
    print(zero_variable_checks)
    # print(stored_check_lookup.shape)
    for check in zero_variable_checks:
        print(stored_check_lookup[check*6:(check+1)*6])
    # print(stored_check_lookup[zero_variable_checks*6])
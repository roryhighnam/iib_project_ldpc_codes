import numpy as np
from galois import GF2
from itertools import permutations

code_directory = '/home/rory/Documents/iib_project_ldpc_codes/hpc_parity_checks/'
code = 4
with open(code_directory+'code_no_'+str(code)+'_n_512_dv_3_dc_6.npy', 'rb') as f:
    parity_check = np.load(f)
parity_check = np.reshape(parity_check, (256,512))

# Check for codewords with minimum distance 2
codewords = []
one_positions = [0,1]
while one_positions[0]<512:
    if one_positions[1] == 512:
        one_positions[0] += 1
        one_positions[1] = one_positions[0]+1
    if one_positions[0] == 511:
        break
    test_codeword = np.zeros(512)
    test_codeword[one_positions] = 1
    codeword = True
    for row in parity_check:
        if np.dot(row, test_codeword) % 2 != 0:
            codeword = False
            break
    if codeword:
        codewords.append(test_codeword)
    one_positions[1] += 1

print(len(codewords))

# Check for codewords with minimum distance 3
codewords = []
one_positions = [0,1,2]
while one_positions[0]<512:

    one_positions[2] += 1
    if one_positions[2] == 512:
        one_positions[1] += 1
        one_positions[2] = one_positions[1] + 1
    


    if one_positions[2] == 512:
        if one_positions[1] != 511:
            one_positions[1] += 1
            one_positions[2] = one_positions[1] + 1
    if one_positions[1] == 512:
        one_positions[0] += 1
        one_positions[1] = one_positions[0] + 1
        one_positions[2] = one_positions[1] + 1
    if one_positions[0] == 511:
        break

    print(one_positions)
    one_positions[2] += 1

#     test_codeword = np.zeros(512)
#     test_codeword[one_positions] = 1
#     codeword = True
#     for row in parity_check:
#         if np.dot(row, test_codeword) % 2 != 0:
#             codeword = False
#             break
#     if codeword:
#         codewords.append(test_codeword)
#     one_positions[1] += 1

# print(len(codewords))
import numpy as np

class regular_LDPC_code():
    '''This class assumes the parity check matrix is in systematic form'''

    def __init__(self,parity_check):
        self.n = len(parity_check[0])
        self.k =  len(parity_check[0]) - len(parity_check)
        self.parity_check = parity_check
        self.generator = np.concatenate((np.eye(self.k), np.transpose(parity_check[:, self.n-self.k:self.n])), axis=1)
        self.dv = np.sum(parity_check, axis=0)[0]
        self.dc = np.sum(parity_check, axis=1)[0]
        self.rate = self.k / self.n

    def encode(self,input_sequence):
        assert len(input_sequence) == self.n, "The input sequence must have the same length as the codewords"


# Example taken from lecture notes
parity_check = np.array([[0,0,1,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,1,0],[0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0],[1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,1,0],
                        [0,1,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1],[1,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0],[0,0,1,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1],
                        [0,0,0,0,1,1,0,1,0,1,0,1,0,0,1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,0],[0,1,0,1,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,1],
                        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0]])

LDPC = regular_LDPC_code(parity_check)
print(LDPC.generator)
print(LDPC.n)
print(LDPC.k)
print(LDPC.dv)
print(LDPC.dc)
print(LDPC.rate)
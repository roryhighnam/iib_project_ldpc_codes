import random
import numpy as np

class BEC:
    '''Class to represent a Binary Erasure Channel. Bit mappings: 0 -> -1, ? -> 0, 1 -> 1. Binary sequences are stored and returned as lists of binary digits.'''

    def __init__(self, erasure_prob):
        self.erasure_prob = erasure_prob
    
    # def transmit(self, input_binary):
    #     output_binary = []
    #     for bit in input_binary:
    #         if random.random() < self.erasure_prob:
    #             output_binary += [0]
    #         else:
    #             output_binary += [-1] if bit == 0 else [1]
    #     return output_binary
    
    def transmit(self, input_binary):
        input_binary = np.where(input_binary==0, -1, input_binary)
        random_vector = np.random.rand(len(input_binary))
        return np.where(random_vector<self.erasure_prob, 0, input_binary)
    
    def new_transmit(self, input_binary):
        random_vector = np.random.rand(len(input_binary))
        return np.where(random_vector<self.erasure_prob, 2, input_binary)
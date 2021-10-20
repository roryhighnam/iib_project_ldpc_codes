import random

class BEC:
    '''Class to represent a Binary Erasure Channel. Bit mappings: 0 -> -1, ? -> 0, 1 -> 1. Binary sequences are stored and returned as lists of binary digits.'''

    def __init__(self, erasure_prob):
        self.erasure_prob = erasure_prob
    
    def transmit(self, input_binary):
        output_binary = []
        for bit in input_binary:
            if random.random() < self.erasure_prob:
                output_binary += [0]
            else:
                output_binary += [-1] if bit == 0 else [1]
        return output_binary
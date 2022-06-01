import numpy as np
import matplotlib.pyplot as plt

def calculate_crit_point(erasure_prob, dv,dc):
    prev_x = 0
    x = 1

    while abs(x-prev_x) > 1e-8:
        prev_x = x
        x = (1/dv) * (dv - 1 - (dv-1)*(1-erasure_prob*x**(dv-1))**(dc-1) + erasure_prob * (dv-1) * (dc-1) * x **(dv-1) * (1 - erasure_prob * x**(dv-1))**(dc-2))
    return x

def calculate_crit_epsilon(dv,dc, low=0.2, high=1):

    while high-low > 1e-8:
        midpoint = 0.5*(high + low)
        if abs(calculate_crit_point(midpoint, 3,6)) < 1e-9:
            low = midpoint

        else:
            high = midpoint
    
    return high

# print(calculate_crit_point(0.3747712850570679, 3, 6))
print(calculate_crit_epsilon(3,6))
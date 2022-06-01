from cgi import test
import csv
import numpy as np
import matplotlib.pyplot as plt
from test_de_threshold import calc_threshold
import scipy.stats as stats
import functools

@functools.cache
def calculate_threshold_y(threshold_erasure,dv,dc):
    previous_value = 0
    test_value = 1
    while np.abs(test_value - previous_value) > 1e-6:
        previous_value = test_value
        test_value = 1 - (1 - threshold_erasure * test_value**(dv-1))**(dc-1)
    return test_value

def calculate_alpha(threshold, dv, dc):
    y_crit = calculate_threshold_y(threshold, dv,dc)
    x_crit = threshold * (y_crit)**(dv-1)
    return threshold * np.sqrt( ((dv-1)/dv) * (1/x_crit - 1/y_crit)) 

# def prob_of_error(n,dv,dc,erasure_prob):
#     threshold = calc_threshold(dv,dc)
#     z = np.sqrt(n) * (threshold - erasure_prob)
#     y = calculate_threshold_y(threshold, dv,dc)
#     x = threshold * (y)**(dv-1)
#     alpha = threshold * np.sqrt( ((dv-1)/dv) * (1/x - 1/y))
#     return stats.norm.cdf(-z/alpha)

erasure_probs = np.linspace(0.32, 0.5, 100)
ns = [5000]
# n = 2000
dv = 3
dc = 6
if __name__ == '__main__':
    threshold = calc_threshold(dv,dc)
    alpha = calculate_alpha(threshold, dv,dc)
    for n in ns:
        # z = np.sqrt(n) * (threshold - erasure_probs - 0.616045 * n**(-2/3))
        z = np.sqrt(n) * (threshold - erasure_probs)
        results = stats.norm.cdf(-z/alpha)
        plt.plot(erasure_probs, stats.norm.cdf(-z/alpha))
    plt.yscale('log')
    plt.show()
n = 5000
# with open('./finite_length_scaling_report/error_prob_n='+str(n)+'.csv', 'w') as csv_file:
#     csv_writer = csv.writer(csv_file)
#     for i in range(len(results)):
#         csv_writer.writerow([erasure_probs[i], results[i]])

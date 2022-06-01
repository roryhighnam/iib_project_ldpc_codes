import matplotlib.pyplot as plt
import numpy as np
import math
import csv
import scipy.linalg as linalg
from test_critical_point_calculator import calculate_crit_point

def f(erasure_prob,dv,dc,y):
    y = (1 - y/(erasure_prob*n))**(1/dv)
    return erasure_prob * y**(dv-1) * (y - 1 + (1-erasure_prob*y**(dv-1))**(dc-1))

def g(erasure_prob,dv,dc,y):
    y = np.array(y)
    return erasure_prob * y**(dv-1) * (y - 1 + (1-erasure_prob*y**(dv-1))**(dc-1)) 

results = []

def check_deg_decoding_path(erasure_prob,dv,dc,degree,y):
    y = np.array(y)
    return math.comb(dc,degree) * (erasure_prob * y**(dv-1))**degree * (1-erasure_prob*y**(dv-1))**(dc-1)

def alt_check_deg_decoding_path(erasure_prob, dv,dc,degree,y):
    y = np.array(y)
    sum = np.zeros_like(y)
    for j in range(degree,dc+1):
        sum += (-1)**(j+degree) * math.comb(j-1, degree-1) * math.comb(dc-1, j-1) * (erasure_prob * y**(dv-2))**j
    return sum

dv = 3
dc = 6
n = 1000
erasure_probs = [0.36, 0.38, 0.4, 0.42, 0.44]
# ys = np.linspace(0,1,1000)

results = []
xs = []
for erasure_prob in erasure_probs:
    ys = np.linspace(0,n*erasure_prob,1000)
    xs.append(ys/(n*erasure_prob))
    result = f(erasure_prob,dv,dc,ys)
    results.append(result[::-1])
    plt.plot(ys/(n*erasure_prob), result[::-1])
    print(erasure_prob)
    minimum_point = n * erasure_prob * (1-(calculate_crit_point(erasure_prob, dv, dc))**dv)
    print((n*erasure_prob - minimum_point)/(n*erasure_prob))
    print("Value at minimum", f(erasure_prob, dv, dc, minimum_point))
    # plt.hlines(dv*n*f(erasure_prob, dv, dc, minimum_point),0,400)
    # plt.vlines(n*erasure_prob - minimum_point, min(0, min(result)), 120)
    plt.plot((n*erasure_prob - minimum_point)/(n*erasure_prob), f(erasure_prob, dv, dc, minimum_point),marker='x' )


plt.show()
# print(results)
# print(xs)

for index,result in enumerate(results):
    with open('./finite_length_scaling_report/peeling_decoder_paths_'+str(erasure_probs[index])+'.csv', 'w') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['x', 'fraction of edges'])

        for i in range(len(result)): 
            csv_writer.writerow([xs[index][i], result[i]])



# erasure_prob = 0.4294
# n = 1000
# ys = np.linspace(0,1,1000)
# results.append(g(erasure_prob, 3,6,ys))
# plt.plot(g(erasure_prob, 3,6,ys))
# for degree in [2,3,4,5,6]:
#     results.append(alt_check_deg_decoding_path(erasure_prob, 3,6,degree,ys)/degree)
#     plt.plot(alt_check_deg_decoding_path(erasure_prob, 3,6,degree,ys)/degree, label='degree = ' + str(degree))

# plt.plot(0.415*2048*ys,2048*3*results[0])
# plt.legend()
# plt.show()

# for point in range(len(ys)):
#     e = 0
#     for degree in range(len(results)):
#         e += results[degree][point]
#     results[0][point] *= results[0][point]/e

# plt.plot(0.415/0.2075*10000*ys, 3*10000*results[0][::-1])

# with open('./finite_length_scaling_report/threshold_results.csv', 'w') as csv_file:
#     csv_writer = csv.writer(csv_file)
#     for index,y in enumerate(ys):
#         data = [y]
#         for result in results:
#             data.append(result[index])
#         csv_writer.writerow(data)


def calculate_initial_values(erasure_prob):
    X0 = np.array([0,0,0,0,0,1])
    print(X0)
    A = np.array([[-1,1,0,0,0,0],[0,-2,2,0,0,0],[0,0,-3,3,0,0],[0,0,0,-4,4,0],[0,0,0,0,-5,5],[0,0,0,0,0,-6]])
    return np.dot(linalg.expm(-A * np.log(erasure_prob)),X0),A

def calculate_distribution(initial_value, dv, A, erasure_prob, time):
    B = np.zeros((6,6))
    B[0,:] = -np.ones(6)
    C = (1/dv)*B + ((dv-1)/dv)*A
    return np.dot(linalg.expm(-C * np.log((1-time)/erasure_prob)), initial_value)


initial_value,A = calculate_initial_values(0.429)
calculate_distribution(initial_value, 3, A, 0.429, 0.7)

results = np.zeros((6,len(ys)))
for y in range(len(ys)):
    if ys[y] > 1 - 0.429 and ys[y] != 1:
        results[:,y] = calculate_distribution(initial_value, 3, A, 0.429, ys[y])
        print(results[:,y])

plt.plot(results[1,:], label='alternative calculation')


plt.legend()
plt.show()
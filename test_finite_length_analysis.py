from asyncore import file_dispatcher
import numpy as numpy
import matplotlib.pyplot as plt
import math
import sympy
from sympy.abc import x
import pickle
import multiprocessing
import os
import functools

# A_dict = {}
# g_dict = {}
# f_dict = {}

data_directory = './finite_length_analysis_data/'

@functools.cache
def f(dc,sigma,dt,coef):
    # if (sigma,dt) in f_dict:
    #     return f_dict[(sigma,dt)][x**coef]
    a = sympy.expand( ((1+x)**(dc-1)-1)**(sigma) * ((1+x)**(dc)-1-dc*x)**dt).as_coefficients_dict()
    # f_dict[(sigma,dt)] = a
    return a[x**coef]

@functools.cache
def A(v,t,dv,dc,N,s,X):
    # First check if value has already been calculated, if so return that
    # if (v,t,s) in A_dict:
    #     return A_dict[(v,t,s)]
    if v==0:
        return 0
    if s==0:
        value = total_stopping_sets(v,dv,dc,t) *(v/N)**X / (math.factorial(v)*math.factorial(dv)**v)
        # A_dict[(v,t,s)] = value
        return value
    
    total = 0
    
    for ds in range(1, dv+1):
        for sigma in range(0, dv-ds+1):
            if s+sigma-ds < 0:
                continue
            omega = (v-1)*dv-(s+sigma-ds)
            fourth_term = math.comb(s+sigma-ds, sigma)
            fifth_term = dc**ds * (ds/s)
            for dt in range(0, math.floor((dv-ds-sigma)/2)+1):
                if dc*(t-dt-sigma)-omega < 0:
                    continue
                first_term = math.comb(dt+ds, dt)
                sixth_term = A(v-1, t-dt-sigma, dv,dc,N,s+sigma-ds,X)
                seventh_term = math.comb(t+s, dt+ds)
                for tau in range(0, dv-ds-sigma-2*dt+1):
                    second_term = f(dc,sigma,dt,dv-ds-tau)
                    third_term = math.comb(dc*(t-dt-sigma)-omega, tau)
                    total += first_term * second_term * third_term * fourth_term * fifth_term * sixth_term * seventh_term
    # A_dict[(v,t,s)] = total
    return total

def B(v,X,k,dv,dc,N):
    sum = 0
    for t in range(k+1):
        for s in range(k-t+1):
            sum += math.comb(k,t+s) * A(v,t,dv,dc,N,s,X)
    sum *= math.factorial(v) * math.factorial(dv)**v
    return sum

def T(v, k, dv, dc):
    return math.factorial(v*dv) * math.comb(k*dc, v*dv)

def p(epsilon, X, n):
    sum = 0
    for e in range(n+1):
        sum += (math.comb(n,e) * epsilon**e * (1-epsilon)**(n-e) * B(e,X)) / T(e)
    return sum

@functools.cache
def g(dc, t):
    if t == 1 or t == 0:
        a = sympy.expand(((1+x)**dc - 1 - dc*x)**t)
        return a
    else:
        a = sympy.expand(g(dc,t-1)*((1+x)**dc - 1 - dc*x))
        return a

def total_stopping_sets(v, dv, dc, t):
    # value = g(dc,t)[x**(v*dv)] * math.factorial(v*dv)
    value = g(dc,t).coeff(x, v*dv) * math.factorial(v*dv)
    return value
    # return g(dc,t)[x**(v*dv)] * math.factorial(v*dv)

def bit_error(epsilon,X,dv,dc,k,n):
    total = 0
    errors = []
    Bs = []
    Ts = []

    for e in range(n+1):
        print("total: ", total)
        print("e: ",e)
        B_value = B(e,X,k,dv,dc,n)
        T_value = T(e,k,dv,dc)
        Bs.append(B_value)
        Ts.append(T_value)
        value = B_value * epsilon**(e) * (1-epsilon)**(n-e) * math.comb(n,e) / T_value
        total += value
        errors.append(value)

    return total,errors,Bs,Ts

n = 1000
erasure_prob = 0.3
dv = 5
dc = 10
X = 0

if __name__ == '__main__':
    total, errors,Bs,Ts = bit_error(erasure_prob, X, dv, dc, int(0.5*n), n)
    print(total)

plt.plot(errors[1:])
plt.yscale('log')
plt.show()

with open(data_directory + 'errors_n=' + str(n) + '_BEC=' + str(erasure_prob) + '_dv=' + str(dv) + '_dc=' + str(dc) + '_X=' +str(X) + '.pickle.nosync', 'wb') as pickle_file:
    pickle.dump(errors, pickle_file)

with open(data_directory + 'B_n=' + str(n) + '_BEC=' + str(erasure_prob) + '_dv=' + str(dv) + '_dc=' + str(dc) + '_X=' +str(X) +  '.pickle.nosync', 'wb') as pickle_file:
    pickle.dump(Bs, pickle_file)

with open(data_directory + 'T_n=' + str(n) + '_BEC=' + str(erasure_prob) + '_dv=' + str(dv) + '_dc=' + str(dc) + '_X=' +str(X) +  '.pickle.nosync', 'wb') as pickle_file:
    pickle.dump(Ts, pickle_file)

# print(T(0,50, 4, 8))

# with open('./A_data.pickle', 'rb') as pickle_file:
#     A_read = pickle.load(pickle_file)
    # print(A_read)

# with open('./A_data_n=' + str(n) + 'BEC=' + str(erasure_prob) + '.pickle.nosync', 'wb') as pickle_file:
#     pickle.dump(A_dict, pickle_file)


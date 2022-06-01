import matplotlib.pyplot as plt

def finite_length_density_evolution_3_6(erasure_prob, iterations,n,threshold=0):
    beta = 0.616949
    epsilon = erasure_prob + beta * n**(-2/3)
    return modified_density_evolution(epsilon, iterations, 3,6,threshold)


def density_evolution(erasure_prob, iterations, dv, dc, threshold=0):
    results = [erasure_prob]
    for it in range(iterations):
        prob = (1-results[-1])**(dc-1)
        prob = erasure_prob * (1-prob)**(dv-1)
        if prob > threshold:
            results.append(prob)
    return results

def modified_density_evolution(erasure_prob, iterations, dv, dc, threshold=0):
    results = [erasure_prob]
    last_prob = erasure_prob
    for it in range(iterations):
        prob = erasure_prob * (1 - (1-last_prob)**(dc-1))**(dv-1)
        new_prob = erasure_prob * (1 - (1-last_prob)**(dc-1))**(dv)

        if new_prob > threshold:
            results.append(new_prob)
        last_prob = prob
    return results

results = density_evolution(0.2, 10, 3, 6, 1e-9)
modified_results = modified_density_evolution(0.3, 11, 3, 6, 2e-50)

# plt.plot(results, label='standard')
# plt.plot(modified_results, label='modified')
# plt.legend()
# plt.show()

# print(modified_results)
# print(finite_length_density_evolution_3_6(0.4, 15, 100, 1e-10))
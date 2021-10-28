import matplotlib.pyplot as plt

def density_evolution(erasure_prob, iterations, dv, dc, threshold):
    results = [erasure_prob]
    for it in range(iterations):
        prob = (1-results[-1])**(dc-1)
        prob = erasure_prob * (1-prob)**(dv-1)
        if prob > threshold:
            results.append(prob)
    return results

erasure_probs = [0.3, 0.4, 0.41, 0.42, 0.43]

for erasure_prob in erasure_probs:
    plt.plot(density_evolution(erasure_prob, 50, 3,6,1e-2))

plt.yscale('log')
plt.show()
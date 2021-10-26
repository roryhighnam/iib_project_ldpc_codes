def density_evolution(erasure_prob, iterations, dv, dc):
    results = [erasure_prob]
    for it in range(iterations):
        prob = (1-results[-1])**(dc-1)
        prob = erasure_prob * (1-prob)**(dv-1)
        results.append(prob)
    return results
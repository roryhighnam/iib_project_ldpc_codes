import matplotlib.pyplot as plt

def density_evolution(erasure_prob, iterations, dv, dc, threshold):
    results = [erasure_prob]
    for it in range(iterations):
        prob = (1-results[-1])**(dc-1)
        prob = erasure_prob * (1-prob)**(dv-1)
        if prob > threshold:
            results.append(prob)
    return results
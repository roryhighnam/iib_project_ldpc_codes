
iterations = 100000
dv = 4
dc = 8
tolerance = 1e-6

def below_threshold(erasure_prob, iterations, dv, dc):
    prob = 1
    last_prob = erasure_prob
    for it in range(iterations):
        last_prob = prob
        prob = (1-last_prob)**(dc-1)
        prob = erasure_prob * (1-prob)**(dv-1)

    return prob < tolerance

def calc_threshold(dv,dc):
    lower_bound_value = 0
    upper_bound_value = 1

    while upper_bound_value-lower_bound_value > 1e-9:
        test_value = 0.5 * (lower_bound_value + upper_bound_value)
        if below_threshold(test_value, iterations, dv, dc):
            lower_bound_value = test_value
        else:
            upper_bound_value = test_value

    return 0.5 * (upper_bound_value + lower_bound_value)



# print(calc_threshold(4,8))
# print(below_threshold(0.43, iterations, 3, 6))

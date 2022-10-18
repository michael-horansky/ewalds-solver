import numpy as np


def scalar_product(vec, k):
    result = vec.copy()
    for i in range(len(vec)):
        result[i] = vec[i] * k
    return(result)

def inner_product(vec1, vec2):
    res = 0.0
    for i in range(len(vec1)):
        res += vec1[i] * vec2[i]
    return(res)

def cross_product(vec1, vec2):
    result = np.zeros(3)
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1] # x component
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2] # y component
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0] # z component
    return(result)

def normalize(vec):
    return(vec / np.sqrt(inner_product(vec, vec)))

def magnitude(vec):
    return(np.sqrt(inner_product(vec, vec)))

def orthogonal_unit_vector(vec1, vec2):
    return(normalize(cross_product(vec1, vec2)))

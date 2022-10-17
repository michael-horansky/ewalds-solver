
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


# read the crystal lattice vectors

print("Input the three primitive lattice vectors in the format 'x y z' (three floats separated by a space) for each vector. Empty string for a_3 will input a unit vector orthogonal to a_1, a_2")

a = [[], [], []]

for i in range(3):
    while(True):
        try:
            cur_a_raw = input(f"  a_{i+1} = ")
            if(i == 2 and cur_a_raw == ''):
                a[2] = orthogonal_unit_vector(a[0], a[1])
                break
            cur_a_raw_list = cur_a_raw.split(' ')
            a[i] = np.array([float(cur_a_raw_list[0]), float(cur_a_raw_list[1]), float(cur_a_raw_list[2])])
            break
        except IndexError:
            print("  Format the input as 'number[space]number[space]number'")
        except ValueError:
            print("  Make sure each part of the input is either an int or a float (exp form permitted)")

# We'll now calculate the primitive reciprocal lattice vectors
# We'll save them as 'reduced vectors', which means we'll omit the factor of 2pi and apply it later

b = [[], [], []]
print("The primitive reciprocal lattice vectors are as follows:")

b_factor = inner_product(a[0], cross_product(a[1], a[2]))
for i in range(3):
    b[i] = scalar_product(cross_product(a[(i+1)%3],a[(i+2)%3]), 1.0/b_factor)
    print(f"  b_{i + 1} = 2 pi ({b[i][0]}, {b[i][1]}, {b[i][2]})")

print("Input the incident wave-vector reduced by the factor of 2pi in the same format.")
k_in = []
while(True):
    try:
        k_in_raw = input(f"  k_in = 2 . pi . ")
        k_in_raw_list = k_in_raw.split(' ')
        k_in = np.array([float(k_in_raw_list[0]), float(k_in_raw_list[1]), float(k_in_raw_list[2])])
        break
    except IndexError:
        print("  Format the input as 'number[space]number[space]number'")
    except ValueError:
        print("  Make sure each part of the input is either an int or a float (exp form permitted)")

print(f"  k_in = 2 pi ({k_in[0]}, {k_in[1]}, {k_in[2]})")

# do a brute force search on the Ewald's sphere
k_out_list = []

zero_threshold = 1e-06

diameter = int(np.ceil(2.0 * magnitude(k_in)))
for h_i in range(- diameter, diameter + 1):
    for h_j in range(- diameter, diameter + 1):
        for h_k in range( - diameter, diameter + 1):
            if h_i == 0 and h_j == 0 and h_k == 0:
                continue
            # G = h_i b_1 + h_j b_2 + h_k b_3, check the Laue condition for viable values of G
            G = (scalar_product(b[0], h_i) + scalar_product(b[1], h_j) + scalar_product(b[2], h_k))
            if (np.absolute(2.0 * inner_product(k_in, G) + inner_product(G, G)) < zero_threshold):
                cur_k_out = k_in + G
                k_out_list.append(cur_k_out)

print("The possible scattered wave-vectors:")
for k_out in k_out_list:
    print(f"  k_out = 2 pi ({k_out[0]}, {k_out[1]}, {k_out[2]})")



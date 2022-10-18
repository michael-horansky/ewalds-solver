

from ewald_defs import *

# A naive algorithm that just checks a box around the Ewalds sphere
# Can fail - if b vectors are too small, it will fail to check every point in the sphere

def reciprocal_lattice_search_naive(b, k_in, diameter, zero_threshold = 1e-08):
    
    k_out_list = []
    points = []
    k_out_matrix = [[], [], []]
    for h_i in range(- diameter, diameter + 1):
        for h_j in range(- diameter, diameter + 1):
            for h_k in range( - diameter, diameter + 1):
            
                # G = h_i b_1 + h_j b_2 + h_k b_3, check the Laue condition for viable values of G
                G = (scalar_product(b[0], h_i) + scalar_product(b[1], h_j) + scalar_product(b[2], h_k))
           
                G_x = G[0]
                G_y = G[1]
                G_z = G[2]
                # print("G_x",G_x)
                # print("G_y",G_y)
                # plt.plot(G_x,G_y,color='red')
                points.append([G_x,G_y,G_z])
                if h_i == 0 and h_j == 0 and h_k == 0:
                    continue
                if (np.absolute(2.0 * inner_product(k_in, G) + inner_product(G, G)) < zero_threshold):
                    cur_k_out = k_in + G
                    k_out_list.append(cur_k_out)
                    k_out_matrix[0].append(cur_k_out[0])
                    k_out_matrix[1].append(cur_k_out[1])
                    k_out_matrix[2].append(cur_k_out[2])
    return(k_out_list, points, k_out_matrix)

# A BFS search that starts at the 0, 0, 0 rec. lattice point and checks branching G vectors until they escape the sphere

def reciprocal_lattice_search_BFS(b, k_in, diameter, zero_threshold = 1e-08):
    
    # Create a queue
    
    queue = [] # each element is a list [h_i, h_j, h_k] to check

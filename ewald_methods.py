
import numpy as np
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

def reciprocal_lattice_search_BFS(b, k_in, diameter, zero_threshold = 1e-08, display_sphere_coef = 1.1):
    
    k_out_list = []
    points = []
    k_out_matrix = [[], [], []]
    
    radius_squared = inner_product(k_in, k_in)
    
    def find_neighbors(node):
        # To reduce redundancy, notice that |x|+|y|+|z| never has to decrease (cause those nodes were already checked). So for each node we only need to check three to five neighbors
        neighbors = []
        node_signs = np.sign(node)
        if node_signs[0] <= 0:
            neighbors.append([ node[0]-1, node[1]  , node[2]   ])
        if node_signs[0] >= 0:
            neighbors.append([ node[0]+1, node[1]  , node[2]   ])
        if node_signs[1] <= 0:
            neighbors.append([ node[0]  , node[1]-1, node[2]   ])
        if node_signs[1] >= 0:
            neighbors.append([ node[0]  , node[1]+1, node[2]   ])
        if node_signs[2] <= 0:
            neighbors.append([ node[0]  , node[1]  , node[2]-1 ])
        if node_signs[2] >= 0:
            neighbors.append([ node[0]  , node[1]  , node[2]+1 ])
        return(neighbors)
    
    def inside_the_sphere(node):
        G = (scalar_product(b[0], node[0]) + scalar_product(b[1], node[1]) + scalar_product(b[2], node[2]))
        k_out = k_in + G
        if inner_product(k_out, k_out) > radius_squared * display_sphere_coef * display_sphere_coef:
            return(False)
        else:
            return(True)
    
    # Create a queue and a list of marked nodes
    
    #queue = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]] # each element is a list [h_i, h_j, h_k] to check
    queue = [[0, 0, 0]]
    marked = [[0, 0, 0]]
    
    while (len(queue) > 0):
        cur_node = queue.pop(0)
        
        cur_G = (scalar_product(b[0], cur_node[0]) + scalar_product(b[1], cur_node[1]) + scalar_product(b[2], cur_node[2]))
        points.append([cur_G[0], cur_G[1], cur_G[2]])
        if cur_node[0] != 0 or cur_node[1] != 1 or cur_node != 2:
            if (np.absolute(2.0 * inner_product(k_in, cur_G) + inner_product(cur_G, cur_G)) < zero_threshold):
                cur_k_out = k_in + cur_G
                k_out_list.append(cur_k_out)
                k_out_matrix[0].append(cur_k_out[0])
                k_out_matrix[1].append(cur_k_out[1])
                k_out_matrix[2].append(cur_k_out[2])
        
        cur_neighbors = find_neighbors(cur_node)
        # ignore already marked neighbors and neighbors outside of the sphere. Queue in the rest.
        for neighbor in cur_neighbors:
            if neighbor in marked:
                continue
            marked.append(neighbor)
            if not inside_the_sphere(neighbor):
                continue
            queue.append(neighbor)
    return(k_out_list, points, k_out_matrix)
            
                
        
        


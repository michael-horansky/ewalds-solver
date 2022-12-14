## JJE - 18/10/2022

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools

from ewald_methods import *
from matplotlib.widgets import Button


# r = 3
# c = 4
# x = np.linspace(0, c, c+1)
# y = np.linspace(0, r, r+1)

# pts = itertools.product(x, y)
# plt.scatter(*zip(*pts), marker='o', s=30, color='red')

# X, Y = np.meshgrid(x, y)
# # deg = np.arctan(Y**3 - 3*Y-X)
# # QP = plt.quiver(X, Y, np.cos(deg), np.sin(deg))
# plt.grid()
# plt.show()

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


points = []
k_out_matrix = [[], [], []] # k[0] is the set of x-components of valid k_out wavevectors etc

# Optimize the search subspace: this will be a box escribed to the sphere
k_in_mag = magnitude(k_in)
diameter = int(np.ceil(2.0 * k_in_mag))
x_min = -k_in[0] - k_in_mag
x_max = -k_in[0] + k_in_mag
y_min = -k_in[1] - k_in_mag
y_max = -k_in[1] + k_in_mag
z_min = -k_in[2] - k_in_mag
z_max = -k_in[2] + k_in_mag


#k_out_list, points, k_out_matrix = reciprocal_lattice_search_naive(b, k_in, diameter, zero_threshold)
k_out_list, points, k_out_matrix = reciprocal_lattice_search_BFS(b, k_in, diameter, zero_threshold, 1.5)

print("The possible scattered wave-vectors:")
for k_out in k_out_list:
    print(f"  k_out = 2 pi ({k_out[0]}, {k_out[1]}, {k_out[2]})")
points = np.array(points)

ax = plt.axes(projection='3d')

u = np.linspace(0, np.pi, 30)
v = np.linspace(0, 2 * np.pi, 30)

x = k_in_mag*np.outer(np.sin(u), np.sin(v)) - k_in[0]
y = k_in_mag*np.outer(np.sin(u), np.cos(v)) - k_in[1]
z = k_in_mag*np.outer(np.cos(u), np.ones_like(v)) - k_in[2]
ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
# Plot the Ewald's sphere

ewald_sphere_plot = ax.plot_surface(x, y, z,alpha=0.2)

# Plot the reciprocal lattice
color_input = 'red'
points_scatter_plot = ax.scatter(points[:,0],points[:,1],points[:,2], s=5, c=color_input)


# Plot the k_in wavevector

k_in_vecotr_plots = ax.quiver(-k_in[0],-k_in[1],-k_in[2],k_in[0],k_in[1],k_in[2],color='black')

# Plot the legit k_out wavevectors

k_out_vecotr_plots = ax.quiver(-k_in[0],-k_in[1],-k_in[2],k_out_matrix[0],k_out_matrix[1],k_out_matrix[2],color='red', linestyles='dashed')


left_bound = -2
right_bound = 2

ax.axes.set_xlim3d(left=left_bound*k_in_mag - k_in[0], right=right_bound*k_in_mag - k_in[0])
ax.axes.set_ylim3d(bottom=left_bound*k_in_mag - k_in[1], top=right_bound*k_in_mag - k_in[1]) 
ax.axes.set_zlim3d(bottom=left_bound*k_in_mag - k_in[2], top=right_bound*k_in_mag - k_in[2]) 



class ToggleObject:
    ind = 0
    def __init__(self, plot_object):
        self.plot_object = plot_object

    def toggle_object(self, event):
        self.ind += 1
        if self.ind % 2 == 0:
             self.plot_object.set_visible(True)
        else:
            self.plot_object.set_visible(False)


callback_scatter = ToggleObject(plot_object=points_scatter_plot)
ax_scatter = plt.axes([0.81, 0.05, 0.1, 0.075])
b_scatter = Button(ax_scatter, 'Points')
b_scatter.on_clicked(callback_scatter.toggle_object)

callback__out_vector = ToggleObject(plot_object=k_out_vecotr_plots)
ax_out_vector = plt.axes([0.70, 0.05, 0.1, 0.075])
b_out_vector = Button(ax_out_vector, 'k out')
b_out_vector.on_clicked(callback__out_vector.toggle_object)

callback_in_vector = ToggleObject(plot_object=k_in_vecotr_plots)
ax_in_vector = plt.axes([0.59, 0.05, 0.1, 0.075])
b_in_vector = Button(ax_in_vector, 'k in')
b_in_vector.on_clicked(callback_in_vector.toggle_object)

callback_ewald_sphere = ToggleObject(plot_object=ewald_sphere_plot)
ax_ewald_sphere = plt.axes([0.48, 0.05, 0.1, 0.075])
b_ewald_sphere = Button(ax_ewald_sphere, 'sphere')
b_ewald_sphere.on_clicked(callback_ewald_sphere.toggle_object)


plt.show()


from get_pmd_data import *
import numpy as np 
import pyvista as pv
import matplotlib.pyplot as plt

pd = pmd_file('small_ssd.pmd')

x_axes = pd.get_header()['axes0'] 
y_axes = pd.get_header()['axes1']
z_axes = pd.get_header()['axes2']

x_axes = np.asarray(x_axes, dtype=float)
y_axes = np.asarray(y_axes, dtype=float)
z_axes = np.asarray(z_axes, dtype=float)

# conversion to km
x_axes *= 0.001
y_axes *= 0.001
z_axes *= 0.001

Nx = len(x_axes)
Ny = len(y_axes)
Nz = len(z_axes)

print('Atmoshpere size:', Nx, Ny, Nz)

T = np.zeros((Nx, Ny, Nz))     # scalar field
B = np.zeros((Nx, Ny, Nz, 3))  # vector field
V = np.zeros((Nx, Ny, Nz, 3))  # vector field

for x in range(Nx):
    print(f"x = {x+1}/{Nx}")
    for y in range(Ny):
        for z in range(Nz):

            node = pd.get_node(ix=x, iy=y, iz=z)

            # scalar
            T[x, y, z] = node['T [K]']

            # vectors
            B[x, y, z, 0] = node['Bx [G]']
            B[x, y, z, 1] = node['By [G]']
            B[x, y, z, 2] = node['Bz [G]']

            V[x, y, z, 0] = node['vx [cm s-1]']
            V[x, y, z, 1] = node['vy [cm s-1]']
            V[x, y, z, 2] = node['vz [cm s-1]']
            

# example: scalar field
grid = pv.RectilinearGrid(x_axes, y_axes, z_axes)
grid["T"] = T.flatten(order="F")
grid.point_data["B"] = B.reshape(-1, 3, order="F")
grid.point_data["V"] = V.reshape(-1, 3, order="F")
grid.save("atmosphere.vtr")


from get_pmd_data import *
import numpy as np 
import pyvista as pv
import matplotlib.pyplot as plt

pd = pmd_file('AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix_conv_KQ_MC.pmd')
field = 'T [K]'

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

data = np.zeros((Nx, Ny, Nz))

for x in range(Nx):
    for y in range(Ny):
        for z in range(Nz):
            data[x, y, z] = pd.get_node(ix=x, iy=y, iz=z)[field]


print('Rendering...')

grid = pv.RectilinearGrid(
    x_axes,
    y_axes,
    z_axes
)

grid["data"] = data.flatten(order="F")

plotter = pv.Plotter()
plotter.add_volume(
    grid,
    scalars="data",
    cmap="inferno",
    opacity="sigmoid"
)
plotter.show()

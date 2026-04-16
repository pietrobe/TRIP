from get_pmd_data import *
import numpy as np 
import pyvista as pv
import matplotlib.pyplot as plt

pd = pmd_file('AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix_conv_KQ_MC.pmd')
field_x = 'vx [cm s-1]'
field_y = 'vy [cm s-1]'
field_z = 'vz [cm s-1]'

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

print('Reading atmoshpere with size:', Nx, Ny, Nz)

V = np.zeros((Nx, Ny, Nz, 3))  # vector field

for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(Nz):
            node = pd.get_node(ix=ix, iy=iy, iz=iz)
            V[ix, iy, iz, 0] = node[field_x]
            V[ix, iy, iz, 1] = node[field_y]
            V[ix, iy, iz, 2] = node[field_z]

print('Rendering...')

grid = pv.RectilinearGrid(x_axes, y_axes, z_axes)

# PyVista wants (Npoints, 3)
vectors = V.reshape(-1, 3, order="F")

grid["V"] = vectors
grid.set_active_vectors("V")


streamlines = grid.streamlines(
    vectors="V",
    source_radius=200000,
    n_points=300
)

plotter = pv.Plotter()
plotter.add_mesh(streamlines, line_width=2, cmap="plasma")

plotter.show_bounds(
    grid='back',
    location='outer',
    all_edges=True
)

plotter.add_axes()
plotter.show()



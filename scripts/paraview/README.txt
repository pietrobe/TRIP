
Execute pmd_to_paraview.py using the right .pmd path in the script. 
Load parview_state.pvsm to inspect the .pmd model in Paraview. 


# Inspect .pmd data:

$ python3
Python 3.9.5 | packaged by conda-forge | (default, Jun 19 2021, 00:32:32) 
[GCC 9.3.0] on linux
>>> from get_pmd_data import *
>>> pd = pmd_file('path_full_or_relative/cai_0B_1V_50x50x133_iter150.pmd')

    Then you can ask for the full header or a single field:

>>> pd.get_header()
{'endian': 0, 'isize': 4, 'fsize': 8, 'version': (2,), 'date': (122, 9, 28, 8, 17, 22), 'period': (1, 1),    ...   'atom_mass': 40.078, 'Aul': 218000000.0, 'E': 4.699722538e-12, 'Jl': 0, 'Ju': 1, 'gl': 0.0, 'gu': 1.0, 'Tref': 6000.0}
>>> pd.get_header('Aul')
218000000.0

    If you want to see just the flags you can make a list of get_header() or input a known wrong field, e.g.:

>>> list(pd.get_header())
['endian', 'isize', 'fsize', 'version', 'date', 'period', 'domain', 'origin', 'dimensions', 'axes0', 'axes1', 'axes2', 'quadrature', 'module', 'comments', 'msize', 'gsize', 'mversion', 'atom_mass', 'Aul', 'E', 'Jl', 'Ju', 'gl', 'gu', 'Tref']

       Axes0,1, and 2 are X, Y, and Z, respectively.

      Regarding the data, I implemented two things. You can read individual nodes:

>>> pd.get_node(ix=7,iy=33,iz=97)
{'epsilon': 6.320206487954581e-05, 'T [K]': 5999.31689453125, 'N [cm-1]': 35.70768038671875,    ...    'continuum opacity [cm-1]': 2.9324800000000064e-16, 'continuum emissivity [cgs]': 5.3167593487302336e-21}

      Or the whole cube (it's only 100 something MB, so it is fine for this one). In this case a legend will be provided to identify the variables:

>>> data, legend = pd.get_cube()
>>> print(data.shape)
(133, 50, 50, 42)

      With the order (nz,ny,nx,variables)

      If we look at the variable in index 1 in the node (ix,iy,iz) = (7,33,97)

>>> print(f'{legend[1]}: {data[97,33,7,1]}')
T [K]: 5999.31689453125

      which you can see that is what was given by pd.get_node() for that coordinate.

      
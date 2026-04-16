# -*- coding: utf-8 -*-

import os,sys,struct
import numpy as np

# Tanaus\'u del Pino Alem\'an
#
# Wrote this quickly to send to Pietro, so it is probably on the top
# of the non-elegant python codes list
#
# A class to manage a (twolevel) pmd file
#

class pmd_file():
    ''' Class to read a pmd file
    '''

    def __init__(self,filename):
        ''' Initialize the class from the system path
        '''

        print('File opening...')

        # Open the file
        try:
            self.__f = open(filename,'rb')
        except FileNotFoundError:
            sys.exit('Could not find the file')
        except:
            raise

        print('File opened')

        valid = self.__get_header()
        if not valid:
            sys.exit('Problem reading header')

        valid = self.__get_mheader() 
        if not valid:
            sys.exit('Problem reading module header')

    def __get_header(self):
        ''' Read the PORTA general header
        '''

        # Try reading header
        try:

            # Initialize size
            self.__hsize = 0

            # Magic label
            label = self.__f.read(8).decode('utf-8')
            if label != 'portapmd':
                print('The file is not a pmd file')
                return False
            self.__hsize += 8

            # Initialize header dictionary
            self.__header = {}

            # Endian
            self.__header['endian'] = struct.unpack('b',self.__f.read(1))[0]
            self.__hsize += 1

            # Int size
            self.__header['isize'] = struct.unpack('b',self.__f.read(1))[0]
            self.__hsize += 1

            # Float size
            self.__header['fsize'] = struct.unpack('b',self.__f.read(1))[0]
            self.__hsize += 1

            # Version
            self.__header['version'] = struct.unpack('i',self.__f.read(4))
            self.__hsize += 4

            # Date
            self.__header['date'] = struct.unpack('i'*6,self.__f.read(24))
            self.__hsize += 24

            # Period
            self.__header['period'] = struct.unpack('bb',self.__f.read(2))
            self.__hsize += 2

            # Domain
            self.__header['domain'] = struct.unpack('ddd',self.__f.read(24))
            self.__hsize += 24

            # Origin
            self.__header['origin'] = struct.unpack('ddd',self.__f.read(24))
            self.__hsize += 24

            # Nodes
            self.__header['dimensions'] = struct.unpack('iii',self.__f.read(12))
            self.__hsize += 12

            # Axes
            for i,n in enumerate(self.__header['dimensions']):
                self.__header[f'axes{i}'] = struct.unpack('d'*8192,self.__f.read(8192*8))
                self.__header[f'axes{i}'] = self.__header[f'axes{i}'][:n]
                self.__hsize += 8192*8

            # Angles
            self.__header['quadrature'] = struct.unpack('ii',self.__f.read(8))
            self.__hsize += 8

            # Module
            self.__header['module'] = self.__f.read(1023).decode('utf-8')
            i = self.__header['module'].find('\x00')
            if i >= 0:
                self.__header['module'] = self.__header['module'][:i]
            self.__hsize += 1023

            # Comments
            self.__header['comments'] = self.__f.read(4096).decode('utf-8').strip()
            i = self.__header['comments'].find('\x00')
            if i >= 0:
                self.__header['comments'] = self.__header['comments'][:i]
            self.__hsize += 4096

            # Module size
            self.__header['msize'] = struct.unpack('i',self.__f.read(4))[0]
            self.__hsize += 4

            # Grid size
            self.__header['gsize'] = struct.unpack('i',self.__f.read(4))[0]
            self.__hsize += 4

        # Could not read
        except:
            raise
            return False

        # Success
        return True

    def __get_mheader(self):
        ''' Read the PORTA twolevel module header
        '''

        if self.__header['module'] != 'twolevel':
            sys.exit('This class is only prepared to read the twolevel module')

        # Try reading header
        try:

            # Initialize size
            self.__mhsize = 0

            # Module version
            self.__header['mversion'] = struct.unpack('i',self.__f.read(4))[0]
            self.__mhsize += 4

            # Atomic mass
            self.__header['atom_mass'] = struct.unpack('d',self.__f.read(8))[0]
            self.__mhsize += 8

            # Einstein coefficient
            self.__header['Aul'] = struct.unpack('d',self.__f.read(8))[0]
            self.__mhsize += 8

            # Transition energy
            self.__header['E'] = struct.unpack('d',self.__f.read(8))[0]
            self.__mhsize += 8

            # Angular momentum
            self.__header['Jl'] = struct.unpack('i',self.__f.read(4))[0]//2
            self.__header['Ju'] = struct.unpack('i',self.__f.read(4))[0]//2
            self.__mhsize += 8

            # Landé factors
            self.__header['gl'] = struct.unpack('d',self.__f.read(8))[0]
            self.__header['gu'] = struct.unpack('d',self.__f.read(8))[0]
            self.__mhsize += 16

            # Doppler width temperature
            self.__header['Tref'] = struct.unpack('d',self.__f.read(8))[0]
            self.__mhsize += 8

            # Bottom boundary temperature
            # Not bothering
            self.__mhsize += 8*self.__header['dimensions'][0]* \
                               self.__header['dimensions'][1]

            # Sanity check
            if self.__mhsize != self.__header['msize']:
                print('Module size not the same than expected in header')
                return False

            # Jump to data
            self.__hjump = self.__hsize + self.__mhsize

        # Could not read
        except:
            raise
            return False

        # Define variable names
        self.__vars = ['epsilon','T [K]','N [cm-1]', \
                       'Bx [G]','By [G]','Bz [G]', \
                       'vx [cm s-1]','vy [cm s-1]','vz [cm s-1]', \
                       'rhoKQ(l)','rhoKQ(u)','JKQ','damping', \
                       'depolarization rate','continuum opacity [cm-1]', \
                       'continuum emissivity [cgs]']
        self.__i0 = 0
        self.__i1 = 8
        self.__rhol = 9
        self.__rhou = 10
        self.__J = 11
        self.__j0 = 12
        self.__j1 = 15

        # Success
        return True

    def get_header(self,flag=None):
        ''' Return the whole header or a specific field
        '''
        if flag is None:
            return self.__header
        else:
            try:
                return self.__header[flag]
            except KeyError:
                print(f'Invalid keyword {flag}')
                print('Valid keywords:',list(self.__header))
            except:
                raise

    def get_node(self,ix,iy,iz):
        ''' Return an array with all variables in a node
        '''
        # Initial jump
        self.__f.seek(self.__hjump,0)

        # Initialize node
        node = {}

        # Compute jump
        jump = self.__header['dimensions'][0]*self.__header['dimensions'][1]* \
               iz*self.__header['gsize']
        jump += self.__header['dimensions'][0]*iy*self.__header['gsize']
        jump += ix*self.__header['gsize']

        # Jump to node
        self.__f.seek(jump,1)

        # Read simple variables
        for i in range(self.__i0,self.__i1+1):
            node[self.__vars[i]] = struct.unpack('d',self.__f.read(8))[0]

        # Get density matrix (lower level)
        node[self.__vars[self.__rhol]] = {}
        node[self.__vars[self.__rhol]]['00'] = struct.unpack('d',self.__f.read(8))[0]
        self.__f.seek(8,1)

        # Get density matrix (upper level)
        node[self.__vars[self.__rhou]] = {}
        node[self.__vars[self.__rhou]]['00'] = struct.unpack('d',self.__f.read(8))[0]
        self.__f.seek(24,1)
        node[self.__vars[self.__rhou]]['10'] = struct.unpack('d',self.__f.read(8))[0]
        self.__f.seek(8,1)
        node[self.__vars[self.__rhou]]['11R'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__rhou]]['11I'] = struct.unpack('d',self.__f.read(8))[0]
        self.__f.seek(32,1)
        node[self.__vars[self.__rhou]]['20'] = struct.unpack('d',self.__f.read(8))[0]
        self.__f.seek(8,1)
        node[self.__vars[self.__rhou]]['21R'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__rhou]]['21I'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__rhou]]['22R'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__rhou]]['22I'] = struct.unpack('d',self.__f.read(8))[0]

        # Get radiation field tensor
        node[self.__vars[self.__J]] = {}
        node[self.__vars[self.__J]]['00'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__J]]['10'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__J]]['11R'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__J]]['11I'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__J]]['20'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__J]]['21R'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__J]]['21I'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__J]]['22R'] = struct.unpack('d',self.__f.read(8))[0]
        node[self.__vars[self.__J]]['22I'] = struct.unpack('d',self.__f.read(8))[0]

        # Read simple variables
        for i in range(self.__j0,self.__j1+1):
            node[self.__vars[i]] = struct.unpack('d',self.__f.read(8))[0]

        # Return
        return node


    def get_cube(self):
        ''' Return an array with all variables in the model
        '''
        # Initial jump
        self.__f.seek(self.__hjump,0)

        # Compute dimension size
        size = self.__header['dimensions'][0]*self.__header['dimensions'][1]* \
               self.__header['dimensions'][2]* \
               (13 + 9 + 2*\
                ((2*self.__header['Jl']+1)**2 + \
                 (2*self.__header['Ju']+1)**2))

        # Get whole data
        data = np.array(struct.unpack('d'*size,self.__f.read(size*8))). \
               reshape((self.__header['dimensions'][2], \
                        self.__header['dimensions'][1], \
                        self.__header['dimensions'][0], \
                        (13 + 9 + 2* \
                         ((2*self.__header['Jl']+1)**2 + \
                          (2*self.__header['Ju']+1)**2))))

        # Legend
        legend = self.__vars[self.__i0:self.__i1+1]
        for K in range(0,2*self.__header['Jl']):
            for Q in range(-K,K):
                legend += [f'Re[rho{K}{Q}(l)]']
                legend += [f'Im[rho{K}{Q}(l)]']
        for K in range(0,2*self.__header['Ju']):
            for Q in range(-K,K):
                legend += [f'Re[rho{K}{Q}(u)]']
                legend += [f'Im[rho{K}{Q}(u)]']
        legend += ['Re[J00]']
        legend += ['Re[J10]']
        legend += ['Re[J11]']
        legend += ['Im[J11]']
        legend += ['Re[J20]']
        legend += ['Re[J21]']
        legend += ['Im[J21]']
        legend += ['Re[J22]']
        legend += ['Im[J22]']
        legend += self.__vars[self.__j0:self.__j1+1]

        # Return
        return data, legend

"""Class for reading grid and input

Reads grid definition from grid_def.all
Reads vertical grid from grid.h5
Reads parameters from input.toml

Jamie Hilditch December 2022
"""

import os
import tomllib

import h5py
import numpy as np
import re

def _parse_grid_def(filepath: str):
    """Parse the grid_def.all file"""

    with open(filepath,'r') as f:
        lines = f.readlines()

    ma = re.search("N[Xx][^0-9]*([0-9]*)",lines[1])
    if ma:
        NX = int(ma.group(1))
    else:
        raise ValueError("Unable to read NX from grid_def.all")
    
    ma = re.search("N[Yy][^0-9]*([0-9]*)",lines[2])
    if ma:
        NY = int(ma.group(1))
    else:
        raise ValueError("Unable to read NY from grid_def.all")

    ma = re.search("N[Zz][^0-9]*([0-9]*)",lines[3])
    if ma:
        NZ = int(ma.group(1))
    else:
        raise ValueError("Unable to read NZ from grid_def.all")

    ma = re.search("N_?[Tt][Hh][^0-9]*([0-9]*)",lines[4])
    if ma:
        NTH = int(ma.group(1))
    else:
        raise ValueError("Unable to read NTH from grid_def.all")
    
    return NX,NY,NZ,NTH
    
def _read_grid_h5(filepath: str):
    """read grid.h5"""
    with h5py.File(filepath,'r') as f:
        G = f['grids/y'][()]
    GF = (G[1:] + G[:-1])/2
    return G,GF 

def _read_input(filepath: str):
    """read input.toml"""
    with open(filepath,'rb') as fp:
        inputs = tomllib.load(fp)
    return inputs

class setup:

    def __init__(self,directory: str = './'):
        self.directory = directory
        self.NX,self.NY,self.NZ,self.NTH = _parse_grid_def(os.path.join(directory,'grid_def.all'))
        self.G,self.GF = _read_grid_h5(os.path.join(directory,'grid.h5'))
        if self.GF.size != self.NY:
            raise ValueError("Size of GF not equal to NY")
        self.inputs = _read_input(os.path.join(directory,'input.toml'))

    @property
    def version(self):
        return self.inputs['VERSION']

    @property
    def scheme(self):
        return self.inputs['SCHEME']

    @property
    def physical(self):
        return self.inputs['PHYSICAL']

    @property
    def timestepping(self):
        return self.inputs['TIMESTEPPING']

    @property
    def output(self):
        return self.inputs['OUTPUT']

    @property
    def initial_conditions(self):
        return self.inputs['INITIAL_CONDITIONS']

    @property
    def forcing(self):
        return self.inputs['FORCING']

    @property
    def velocity_bcs(self):
        return self.inputs['VELOCITY_BCS']

    @property
    def scalars(self):
        return self.inputs['SCALARS']

    def get_scalar(self,index: int):
        return self.inputs['SCALARS'][index]

    @property
    def x_grid(self):
        return np.linspace(0,self.physical['LX'],self.NX,endpoint=False)

    @property
    def z_grid(self):
        return np.linspace(0,self.physical['LZ'],self.NZ,endpoint=False)



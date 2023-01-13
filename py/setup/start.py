"""Class for creating a start.h5 file

Mimic an out.%%%%%%.h5 file but with initial conditions
This allows for more complicated initial conditions than can be 
easily achieved by editing the fortran code

Jamie Hilditch December 2022
"""
import os 

import dask.array as da
import h5py
import numpy as np

from setup import setup

def _interpolate_to_fractional_grid(V):
    return (V[:,:-1,:] + V[:,1:,:])/2

class start():

    def __init__(self,setup: setup):
        self.setup = setup
        self.directory = setup.directory
        self.NX = setup.NX
        self.NY = setup.NY 
        self.NZ = setup.NZ
        self.NTH = setup.NTH
        self.U = da.empty((setup.NZ,setup.NY,setup.NX))
        self.V = da.empty((setup.NZ,setup.NY+1,setup.NX))
        self.W = da.empty((setup.NZ,setup.NY,setup.NX))
        self.TH = da.empty((self.NTH,setup.NZ,setup.NY,setup.NX))
        self.x = setup.x_grid 
        self.G = setup.G 
        self.GF = setup.GF 
        self.z = setup.z_grid
        self.ZF, self.YF, self.XF = np.meshgrid(setup.z_grid,setup.GF,setup.x_grid,sparse='true',indexing='ij')
        self.ZB, self.YB, self.XB = np.meshgrid(setup.z_grid,setup.G,setup.x_grid,sparse='true',indexing='ij')
        self.Time = 0
        self.Save_Flow_Time = setup.output['SAVE_FLOW_DT']
        self.Save_Stats_Time = setup.output['SAVE_STATS_DT']
        self.Save_Movie_Time = setup.output['SAVE_MOVIE_DT']
        self.Time_Step = 0

    def write_start_file(self):
        filepath = os.path.join(self.directory,'start.h5')
        with h5py.File(filepath,'w') as f:
            f.attrs.create('Resolution',[self.NX,self.NY,self.NZ],dtype=np.dtype('i'))
            grp = f.create_group('Timestep')
            grp.attrs.create('Time',self.Time,dtype=np.float64)
            grp.attrs.create('Save_Flow_Time',self.Save_Flow_Time,dtype=np.float64)
            grp.attrs.create('Save_Stats_Time',self.Save_Stats_Time,dtype=np.float64)
            grp.attrs.create('Save_Movie_Time',self.Save_Movie_Time,dtype=np.float64)
            grp.attrs.create('Time_Step',self.Time_Step,dtype=np.dtype('i'))

        hdf5_variables = {
            '/Timestep/U': self.U,
            '/Timestep/V': _interpolate_to_fractional_grid(self.V),
            '/Timestep/W': self.W
        }
        for i in range(self.NTH):
            hdf5_variables[f'/Timestep/TH{i+1}'] = self.TH[i,:,:,:]

        da.to_hdf5(filepath,hdf5_variables)
            
if __name__ == "__main__":
    se = setup('./2d_unforced')
    star = start(se)
    star.write_start_file()
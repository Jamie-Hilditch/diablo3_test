"""Grid generation classes"""
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

# Define an abstract class for generating grids
# Require a description of the grid and a function that turns 
# scaled indices [-1/N,0/N,...,N/N,(N+1)/N] into grid points
# Also implement a default empty parameter information property
class GridGenerator(ABC):

    @property
    @abstractmethod
    def description(self) -> str:
        pass

    @property
    def parameter_info(self) -> str:
        return ""

    @abstractmethod
    def create_grid(self,J: NDArray[np.float64]) -> NDArray[np.float64]:
        pass

# Abstract class for stretched grid generators with a stretching parameter CS
class StretchedGridGenerator(GridGenerator):
    def __init__(self,CS: float):
        self.CS = CS 

    @property
    def parameter_info(self) -> str:
        return f"Stretching Parameter CS = {self.CS}"

# Define the different grid types
class uniform_positive(GridGenerator):
    description = "Uniform grid [0,1]"
    def create_grid(self,J: NDArray[np.float64]) -> NDArray[np.float64]:
        return J - 1
    
class uniform_symmetric(GridGenerator):
    description = "Uniform grid [-1/2,1/2]"
    def create_grid(self, J: NDArray[np.float64]) -> NDArray[np.float64]:
        return J - 1/2

class uniform_negative(GridGenerator):
    description = "Uniform grid [-1,0]"
    def create_grid(self, J: NDArray[np.float64]) -> NDArray[np.float64]:
        return J - 1

class high_resolution_ends_positive(StretchedGridGenerator):
    description = "Tanh stretching for high resolution at ends [0,1]"
    def create_grid(self, J: NDArray[np.float64]) -> NDArray[np.float64]:
        return (np.tanh(self.CS*(2*J-1))/np.tanh(self.CS) + 1)/2

class high_resolution_ends_symmetric(StretchedGridGenerator):
    description="Tanh stretching for high resolution at ends [-1/2,1/2]"
    def create_grid(self, J: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.tanh(self.CS*(2*J-1.0))/np.tanh(self.CS)/2

class high_resolution_centre_symmetric(StretchedGridGenerator):
    description="Cubic stretching for high resolution in the centre [-1/2,1/2]"
    def create_grid(self, J: NDArray[np.float64]) -> NDArray[np.float64]:
        return (self.CS*(2*J-1.0)**3+2*J-1.0)/(self.CS+1)/2

# Enumerate and store all the grid types using a dict 
GRID_GENERATORS: dict[ int, GridGenerator ] = {
    1: uniform_positive,
    2: uniform_symmetric,
    3: uniform_negative,
    4: high_resolution_ends_positive,
    5: high_resolution_ends_symmetric,
    6: high_resolution_centre_symmetric
}
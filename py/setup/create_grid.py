"""Create a grid for diablo simulations

Grid types are based off John Taylor's matlab code (2005/10/23)
Requires h5py, matplotlib, and numpy

J. Hilditch 2022/12/06
"""

import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

import grid_generators as gg

CS_DEFAULT = 1.5
PLOTFILE_DEFAULT = "grid.png"
OUTFILE_DEFAULT = "grid.h5"

def create_grid(N: int, H: float, generator: gg.GridGenerator): 
    """Create a grid 
    Arguments:
        N: The number of grid points - must match grid_def.all
        H: Height of the domain
        generator: Implementation of the GridGenerator class

    Outputs:
        G: the locations of points on the base grid
        GF: the locations of points on the fractional grid

    The grid definitions are produced as follows:
    1 First, define the G grid based on a specified function which
      loops over j=0,N+2 where 1 and N+1 correspond to the bottom
      and top wall values respectively.
    2 Then, define the fractional grid GF halfway between neighboring
      G grid locations.
    3 Stretch both GF and G so that GF(1) and GF(N) now correspond
      to the upper and lower walls.
    """


    # First make a scaled index vector
    J = np.arange(-1,N+2,1,dtype=np.float64)/N

    # create the base grid
    G = H*generator.create_grid(J)

    # Now generate the fractional grid
    GF = (G[:-1] + G[1:])/2

    # Rescale to place GF[1] and GF[N] at the walls
    gf_low = GF[1]
    gf_up = GF[N]
    g_low = G[1]
    g_up = G[N+1]

    GF = GF*(g_up - g_low)/(gf_up-gf_low)
    G = G*(g_up - g_low)/(gf_up-gf_low)

    # Shift grids
    gf_low = GF[1]
    shift = g_low - gf_low

    GF = GF + shift
    G = G + shift

    # Trim G to remove excess point on each end
    G = G[1:-1]

    return G,GF

_GRID_TYPE_HELP_STR = (
    "Available grid types are (default = 1):\n  " + 
    '\n  '.join((f'{i}) {generator.description}' for i,generator in gg.GRID_GENERATORS.items()))
)

_CS_HELP_STR = (
    "Grid stretching parameter for use with grid types " + 
    ', '.join(str(i) for i,generator in gg.GRID_GENERATORS.items() if issubclass(generator,gg.StretchedGridGenerator)) + 
    f" (default = {CS_DEFAULT})"
)

def positive_integer(arg: str):
    x = int(arg)
    if x <= 0:
        raise ValueError("Argument must be strictly positive")
    return x

def construct_generator(args: argparse.Namespace) -> gg.GridGenerator:
    """Constructs a GridGenerator given the input arguments"""

    generator_class = gg.GRID_GENERATORS[args.type]

    # pass appropriate command line arguments to the grid generator class
    if issubclass(generator_class,gg.StretchedGridGenerator):
        generator = generator_class(CS=args.CS)
    elif issubclass(generator_class,gg.GridGenerator):
        generator = generator_class()
    else:
        raise ValueError('Grid generation class must inherit from GridGenerator')
    return generator

def plot_grid(G: NDArray[np.float64],GF: NDArray[np.float64],plotfile: str):
    MS = 1
    fig, ax = plt.subplots(figsize=(12,4),nrows=2,sharex=True,constrained_layout=True)
    ax[0].plot(GF[1:-1],np.diff(G),'ko',ms=MS)
    ax[0].set_ylabel('DY - defined at GF points')
    ax[0].axvline(GF[1],color='grey',zorder=-1)
    ax[0].axvline(GF[-2],color='grey',zorder=-1)
    ax[1].plot(G,-0.5*np.ones_like(G),'ko',ms=1)
    ax[1].plot(GF,0.5*np.ones_like(GF),'ro',ms=1)
    ax[1].set_xlabel('y - Base grid (black), Fractional grid (red)')
    ax[1].yaxis.set_visible(False)
    ax[1].set_ylim([-1,1])
    ax[1].axvline(GF[1],color='grey',zorder=-1)
    ax[1].axvline(GF[-2],color='grey',zorder=-1)
    fig.savefig(plotfile or PLOTFILE_DEFAULT)
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description = 'Create the vertical grid for a diablo simulation',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('NY',type=positive_integer, help="Number of grid points in the fractional grid")
    parser.add_argument('-t','--type',type=int, choices=gg.GRID_GENERATORS.keys(), default=1,
        metavar='Grid Type', help=_GRID_TYPE_HELP_STR)
    parser.add_argument('-H','--height',type=float, default=1, help="Domain height (default = 1)")
    parser.add_argument('-c','--CS',type=float, default=CS_DEFAULT, help=_CS_HELP_STR)
    parser.add_argument('-o','--outfile',type=str, default=OUTFILE_DEFAULT, 
        help='HDF5 file to save grid to (default grid.h5)')
    parser.add_argument('-p','--plot',action='store_true', help='Save a plot of the grid')
    parser.add_argument('--plotfile',type=str,default='',
        help=('File to save plot of the grid to (default grid.png).\n' + 
            'Setting this automatically sets plot to true.'))


    args = parser.parse_args()
    
    generator = construct_generator(args)

    print("Generating a grid of type:\n  " + generator.description)
    print(f"With {args.NY} points and height {args.height}")
    if generator.parameter_info:
        print("Parameter Information:\n  " + generator.parameter_info)

    # create the grid
    G,GF = create_grid(N=args.NY,H=args.height,generator=generator)

    # compute and write out the maximum stretching ratio
    dGF = np.diff(GF)
    GFr = dGF[1:]/dGF[:-1]
    print(f"Maximum grid stretching ratio: {np.max(GFr):.2f}")

    # save the grid
    with h5py.File(args.outfile,'w') as f:
        ygrid = f.create_dataset('grids/y',(args.NY+1,), dtype=np.float64)
        ygrid[:] = G


    # plot the grid
    if args.plot or args.plotfile:
        plot_grid(G,GF,args.plotfile)
        


if __name__ == "__main__":
    main()
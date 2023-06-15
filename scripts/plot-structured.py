import sys
import argparse

ADIOS2=False
HDF5=False
try:
    import adios2
    ADIOS2=True
    print("Loaded ADIOS2")
except ImportError:
    print("Failed to load ADIOS2 interface")
    print("Attempting to continue with HDF5")
    import h5py
    HDF5=True

import numpy as np
import matplotlib
# Change matplotlib backend
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Define argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Input data file name",
        type=str)
parser.add_argument("-p", "--prefix", help="Output prefix",
        type=str)

def plot_field(field, fig, axs, index, nlevels, name):
    cf = axs[index].contourf(field, levels=nlevels)
    axs[index].set_title(name)

    divider = make_axes_locatable(axs[index])
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    fig.add_axes(ax_cb)

    plt.colorbar(cf, cax=ax_cb)
    ax_cb.yaxis.tick_right()

def read_hdf5(fname):
    with h5py.File(fname, "r") as data:
    
        p = data["Step0"]["p"][()]
        u = data["Step0"]["u"][()]
        v = data["Step0"]["v"][()]

    return u, v, p

def read_adios2(fname):
    with adios2.open(fname, "r") as data:
        for fstep in data: # Expecting only one step

            # For some reason variable names are prefixed with "/"
            p = fstep.read("/p")
            u = fstep.read("/u")
            v = fstep.read("/v")

            break

    return u, v, p

def reshape_data(var, cps):
   
   # Reshape data
   VAR = np.reshape(var, (cps, cps))

   return VAR

def get_grid(var):
    
   cps = int(len(var)**(1/2)) # Cells per side (assumed square mesh)
   L   = 1.0 # Length of side (assumed square mesh)
   dx  = L / cps
   x   = np.linspace(0.5 * dx, L - 0.5 * dx, num=cps)

   return x, cps

def main():

    # Load data
    if len(sys.argv) > 3:
        fname = parser.file
        outstub = parser.prefix
    else:
        fname = "LidDrivenCavity.sol.h5"
        outstub = ""

    if (ADIOS2):
        u, v, p = read_adios2(fname)
    elif (HDF5):
        u, v, p = read_hdf5(fname)

    x, cps = get_grid(p)
    P = reshape_data(p, cps)
    U = reshape_data(u, cps)
    V = reshape_data(v, cps)
        
    mp  = int(np.floor((cps + 1) / 2))
    
    # Plot
    nlevels = 100
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    fig.tight_layout(pad=2.5, w_pad=2.5, h_pad=2.5)
    
    plot_field(P, fig, axs, 0, nlevels, "Pressure")
    plot_field(U, fig, axs, 1, nlevels, "u")
    plot_field(V, fig, axs, 2, nlevels, "v")
    
    filename = outstub + f"contours_{cps:04d}.png"
    fig.savefig(filename, dpi = 200)
    plt.close(fig)
    
    # Ghia data (Re=100)
    
    # u - plotted through vertical centreline
    yu = [0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1]
    ug = [0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1]
    
    xv = [0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.96099, 0.9688, 1]
    vg = [0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0]
    
    fig, (ax_u, ax_v) = plt.subplots(1, 2, figsize = (10, 5))
    ax_u.plot(x, U[:,mp], color = 'blue', label = 'ccs')
    ax_u.plot(yu, ug, color = 'red', label = 'ghia')
    ax_u.set_title("u")
    
    ax_v.plot(x, V[mp,:], color = 'blue', label = 'ccs')
    ax_v.plot(xv, vg, color = 'red', label = 'ghia')
    ax_v.set_title("v")
    
    ax_u.legend()
    fig.tight_layout(pad=2.5, w_pad=2.5, h_pad=2.5)
    fig.savefig(outstub + f"slices_{cps:04d}.png", dpi = 200)

if __name__ == "__main__":
    main()

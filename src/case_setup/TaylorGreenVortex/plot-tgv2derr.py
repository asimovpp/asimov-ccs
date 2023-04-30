"""
A script to plot the error log from TGV2D runs.
Across a range of mesh refinements (as defined by M) the final error in each variable is plotted
alongside theoretical O1 and O2 convergence rates.
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

fmt = ScalarFormatter()

M=[
        16,
        32,
        64,
        128,
        256
]

def read_last(filename):
    with open(filename, "r") as f:
        for row in f:
            pass
        lastline = row
    return row

def parse_error(filename):
    
    # Want errors at end of simulation
    line = read_last(filename)
    
    words = line.split()
    u = float(words[1])
    v = float(words[2])

    return u, v

def ideal(h, y0, order=1):

    id = []
    for hi in h:
        id.append(y0 * (hi / h[0])**order)

    return id

def plot_error(h, err, varname):

    plt.plot(h, err, label=varname,
             ls="", marker="o")
    plt.plot(h, ideal(h, err[0], 1), label="O1")
    plt.plot(h, ideal(h, err[0], 2), label="O2")

    plt.xscale("log", base=2)
    plt.yscale("log")

    plt.xlabel("h (relative)")
    plt.ylabel("Error RMS")
    plt.legend()

    ax = plt.gca()
    ax.xaxis.set_major_formatter(fmt)
    
    figname = "err-" + varname + ".pdf"
    plt.savefig(figname, bbox_inches="tight")

    plt.close()

def relative_h(M):

    h = []
    for m in M:
        h.append(M[0] / m)
    return h

def main():
    
    uh = []
    vh = []
    for m in M:
        filename = "tgv2d-err." + str(m) + ".log"
        u, v = parse_error(filename)

        uh.append(u)
        vh.append(v)

    h = relative_h(M)
    plot_error(h, uh, "u")
    plot_error(h, vh, "v")

if __name__ == "__main__":
    main()

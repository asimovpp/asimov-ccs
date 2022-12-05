"""
A script to plot the error log from TGV2D runs.
Across a range of mesh refinements (as defined by M) the final error in each variable is plotted
alongside theoretical O1 and O2 convergence rates.
"""

import matplotlib.pyplot as plt

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

def ideal(x, y0, order=1):

    id = []
    for xh in x:
        id.append(y0 * (x[0] / xh)**order)

    return id

def plot_error(x, err, varname):

    plt.plot(x, err, label=varname)
    plt.plot(x, ideal(x, err[0], 1), label="O1")
    plt.plot(x, ideal(x, err[0], 2), label="O2")

    plt.xscale("log", base=2)
    plt.yscale("log")

    plt.xlabel("N")
    plt.ylabel("Error RMS")
    plt.legend()

    figname = "err-" + varname + ".pdf"
    plt.savefig(figname, bbox_inches="tight")

    plt.close()

def main():
    
    uh = []
    vh = []
    for m in M:
        filename = "tgv2d-err." + str(m) + ".log"
        u, v = parse_error(filename)

        uh.append(u)
        vh.append(v)

    plot_error(M, uh, "u")
    plot_error(M, vh, "v")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import argparse
import itertools
import re
import numpy as np
from matplotlib import pyplot
from matplotlib import rcParams
from matplotlib import rc
from matplotlib.ticker import FuncFormatter

parser = argparse.ArgumentParser(description='Displays residuals file')

parser.add_argument('path_or_jobid', type=str, nargs='?', help='jobid to display the residuals of, or path to residuals.log', default='.')
parser.add_argument('--l2', action='store_true', default=False,
                    help='Plots L2 norm instead of Linf')
parser.add_argument('--both-l', action='store_true', default=False,
                    help='Plots both Linf and L2 norms')
parser.add_argument('--both-l-sidebyside', action='store_true', default=False,
                    help='Plots both Linf and L2 norms, side by side')
parser.add_argument('--max-it',  type=int, default=40000,
                    help='Maximum number of iterations to display')
parser.add_argument('--end', action='store_true', default=False,
                    help='Display only the end of the residuals')
parser.add_argument('--very-end', action='store_true', default=False,
                    help='Display only the last 400 iterations of the residuals')

args = parser.parse_args()
disp_L2 = args.l2
disp_Linf = not args.l2

disp_sidebyside = args.both_l_sidebyside 

if disp_sidebyside or args.both_l:
    disp_L2 = True
    disp_Linf = True

DISPLAY_END = args.end
MAX_line = args.max_it
if args.very_end:
    MAX_line = 400
    DISPLAY_END = True



def get_wdir(jobid):
    """
        Get the working directory from the job id using scontrol (slurm)
    """
    jobid_str = str(jobid)
    proc = subprocess.Popen(["scontrol", "show", "job", jobid_str], stdout=subprocess.PIPE)
    output = proc.communicate()[0]
    if proc.returncode != 0:
        raise ValueError('Error with scontrol')
    workdir = re.findall("WorkDir=([^\s]*)", str(output))[0]
    return workdir



# get file to load
foldername = ''
if args.path_or_jobid.isdigit():
    jobid = int(args.path_or_jobid)
    foldername = get_wdir(jobid)
else:
    if os.path.isfile(args.path_or_jobid):
        filename = args.path_or_jobid
    elif os.path.isdir(args.path_or_jobid):
        foldername = args.path_or_jobid
    else:
        raise ValueError('File not found')

if foldername:
    if os.path.isfile(foldername+"/residuals.log"):
        filename = foldername+"/residuals.log"
    else:
        raise ValueError('File not found')


def load_header(filename):
    """
        Loads residuals file header and compute size and ith_line
            to make sure only useful columns and a limited amount of lines are loaded
    """

    with open(filename) as fin:
        string_header = fin.readline().rstrip()
        first_line = np.fromstring(fin.readline(), sep=" ")
        for i, l in enumerate(fin):
            pass
        file_length = i + 2


    # select only useful columns
    column2extract = [0,1]

    header = re.findall('\"([^\"]+)\"', string_header)
    header = string_header.split()
    for idx,_ in enumerate(header):
        if ((header[idx].startswith("L2_") and disp_L2) or \
                (header[idx].startswith("Linf_") and disp_Linf)):
            column2extract.append(idx)

    header = [ header[idx] for idx in column2extract ]

    # Load either then end of the file or every ith line
    if DISPLAY_END:
        ith_begining = max(0,file_length - MAX_line)
        ith_line = 1
    else:
        ith_begining = 0
        ith_line = max(1,int(1+(file_length / (MAX_line+1))))
        if ith_line > 10:
            print("Warning: you are plotting only each "+str(ith_line)+" iterations")

    return (header, file_length, ith_begining, ith_line, column2extract)



def load_residuals(filename, ith_begining, ith_line, column2extract):
    """
        Load part of the residuals file, from ith_begining to the end every ith_line lines
    """

    with open(filename) as fin:
        residuals_return = np.loadtxt(itertools.islice(fin, ith_begining, None, ith_line), dtype=float, usecols=column2extract, skiprows=1)

    return residuals_return


# Read of the residual.dat
header, file_length, ith_begining, ith_line, column2extract = load_header(filename)
residuals = load_residuals(filename, ith_begining, ith_line, column2extract)
print("Plot residuals for:   " + filename)
print("Data point displayed: " + str(len(residuals)) + " (out of " + str(file_length) + ")")

##################################
# Start display related section
##################################

rcParams['font.family'] = 'serif'
rcParams["axes.grid"] = True

# color palette definition (V2 from  https://matplotlib.org/users/dflt_style_changes.html#colors-color-cycles-and-color-maps)
colorlist = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
             '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
             '#bcbd22', '#17becf']

if disp_sidebyside:
    f, (ax1,ax2) = pyplot.subplots(1, 2, sharey=True)
else:
    f,ax1 = pyplot.subplots()

# if len(filename) >= 50:
#     f.canvas.set_window_title('...' + filename[-50:])
# else:
#     f.canvas.set_window_title(filename)

f.canvas.mpl_connect('close_event', lambda evt: sys.exit(0))

colorid = 0
usedcolor = {}
plot_lines = [[]]*len(header)

# set up legend, line style etc.
for idx, column in enumerate(header):

    if not column.startswith("L"):
        continue

    fieldname = re.sub("^[^_]*_","", column)
    
    linealpha = 1
    if fieldname in usedcolor:
        if not disp_sidebyside:
            linealpha = 0.6
    else:
        usedcolor[fieldname] = colorlist[colorid]
        colorid = ( colorid + 1) % len(colorlist)

    linecolor = usedcolor[fieldname]

    if disp_sidebyside:
        label = re.sub("_", " ", column)
        label = re.sub("L2 ", "", label)
        label = re.sub("Linf ", "", label)
    else:
        label = re.sub("_", " ", column)
        label = re.sub("L2 ", "L_2 - ", label)
        label = re.sub("Linf ", "L_\\\infty - ", label)

    if idx <= (len(header)/2) or not disp_sidebyside:
        ax = ax1
    else:
        ax = ax2

    plot_lines[idx], = ax.semilogy([], [], alpha=linealpha, color=linecolor, label=r"$" + label + "$")


def fill_plot(residuals, header, plot_lines):
    """
        Fill the plot with residuals data
    """
    for idx, column in enumerate(header):

        if not column.startswith("L"):
            continue

        plot_lines[idx].set_xdata(range(len(residuals[:,0])))
        plot_lines[idx].set_ydata(residuals[:,idx])


def draw_plot():
    """
        autoscale axis and draw plot
    """
    ax1.relim()
    ax1.autoscale_view()

    if disp_sidebyside:
        ax2.relim()
        ax2.autoscale_view()

    pyplot.draw()


def format_tick_labels(x, pos):
    #idx = (np.abs(residuals[:,0] - x)).argmin()
    try:
        timestep = int(residuals[int(x),0])
        return "$ {} $\n$({})$".format(int(x), timestep)
    except:
        return "$ {} $".format(int(x))


# Axis label and titles
ax1.xaxis.set_major_formatter(FuncFormatter(format_tick_labels))
ax1.set_xlabel("Iteration (Timestep)")
if disp_sidebyside:
    ax2.xaxis.set_major_formatter(FuncFormatter(format_tick_labels))
    ax2.set_xlabel("Iteration (Timestep)")

ax1.set_ylabel("Residuals")
if disp_sidebyside:
    ax1.set_title("$L_2$")
    ax2.set_title("$L_{\infty}$")

leg = pyplot.legend()
leg.get_frame().set_linewidth(0.0)

fill_plot(residuals, header, plot_lines)
draw_plot()
pyplot.show()

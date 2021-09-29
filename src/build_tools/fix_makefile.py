import sys
import os

infilename = sys.argv[1]
outfilename = sys.argv[2]

with open(infilename, "r") as infile, open(outfilename, "w") as outfile :
    for line in infile:
        outfile.write(line)
        outfile.write("\t$(FC) $(FFLAGS) -o $@ -c $< $(INC)\n")

import sys
import os

infilename = sys.argv[1]
outfilename = sys.argv[2]

with open(infilename, "r") as infile, open(outfilename, "w") as outfile :
    for line in infile:
        inrule = line.split(" ")
        outrule = ["obj/" + os.path.basename(x) if x.endswith(".o") else x for x in inrule]
        outfile.write(" ".join(outrule))
	

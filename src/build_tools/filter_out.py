import sys
import os

ignore = sys.argv[1].split(" ")
files = sys.argv[2].split(" ")
out = [x for x in files if os.path.basename(x) not in ignore]
print(" ".join(out))

import matplotlib.pyplot as plt
import sys

xs = []
ys = []
#input is either csr_orig.txt or csr_new.txt
with open(sys.argv[1], "r") as f:
  i = 0
  for l in f:
    for nb in map(int, filter(None, l.split(" "))):
      ys.append(i)
      xs.append(nb)
    i = i + 1

plt.scatter(xs, ys, marker=",")
plt.show()

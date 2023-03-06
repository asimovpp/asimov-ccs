import re
import sys

def parse_file(filename, data, categories):
  with open(filename, "r") as f:
    for line in f:
      match = re.search("@build (.+)", line)
      if match:
        m = match.group(1).strip().upper()
        if m in categories:
          data[m].append(filename.replace(".f90", ".o"))

if __name__ == "__main__":
  suffix = "_OBJ"
  obj_categories = ["CAF", "MPI", "MPI+OMP", "PETSC", "KERNEL"]

  data = {x:[] for x in obj_categories}
  for filename in sys.argv[1:]:
    parse_file(filename, data, obj_categories)

  for k,v in data.items():
    print(k + suffix, "=", " ".join(set(v)))

  sys.exit(0)

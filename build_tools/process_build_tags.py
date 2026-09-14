import re
import sys
import os

def parse_file(filename, data, categories, obj_dir):
  with open(filename, "r") as f:
    for line in f:
      match = re.search("@build (.+)", line)
      if match:
        m = match.group(1).strip().upper()
        if m in categories:
          filename = filename.replace(".f90", ".o")
          filename = obj_dir + "/" + os.path.basename(filename)
          data[m].append(filename)

if __name__ == "__main__":
  suffix = "_OBJ"
  obj_categories = ["CAF", "MPI", "MPI+OMP", "PETSC"]
  obj_dir = sys.argv[1]

  data = {x:[] for x in obj_categories}
  for filename in sys.argv[2:]:
    parse_file(filename, data, obj_categories, obj_dir)

  for k,v in data.items():
    print(k + suffix, "=", " ".join(set(v)))

  sys.exit(0)

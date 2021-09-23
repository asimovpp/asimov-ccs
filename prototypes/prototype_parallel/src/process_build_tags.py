import re, sys

def parse_file(filename, data, categories):
  with open(filename, "r") as f:
    for line in f:
      match = re.search("@build (.+)", line)
      if (match and match.group(1).strip() in categories):
        data[match.group(1).strip()].append(filename.replace(".f90", ".o"))

if __name__ == "__main__":
  suffix = "_OBJ"
  obj_categories = ["CAF", "MPI", "MPI+OMP"]
  
  data = {x:[] for x in obj_categories}
  for filename in sys.argv[1:]:
    parse_file(filename, data, obj_categories)
 
  for k,v in data.items():
    print(k + suffix, "=", " ".join(set(v)))
  
  sys.exit(0)

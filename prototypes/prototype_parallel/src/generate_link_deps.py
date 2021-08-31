import sys
import os
import json

# define conflicting configurations
conflicts = \
[
]

# define how to map config options to filenames
method_file_mapping = \
{
  "parallel_env": {"mpi": "parallel_env_mpi",
                   "caf": "parallel_env_caf"},
  "parallel_utils": {"mpi": "parallel_utils_mpi",
                     "caf": "parallel_utils_caf"},
  "compute": {"mpi":     "compute_mpi",
              "mpi_omp": "compute_mpi_omp",
              "caf":     "compute_caf"},
  "collectives": {"mpi": "parallel_collectives_mpi",
                  "none": None},
  "parallel_types": {"mpi": "parallel_types_mpi",
                     "none": None},
}


if __name__ == "__main__":
  with open(sys.argv[1]) as f:
    config = json.load(f)
  print("config read: ", config)

  print("checking config for conflicts")
  for c in conflicts:
    # check for conflicts that appear as subsets: of config
    if (c.items() <= config.items()):
      print("config has following illegal combination of settings: ", c)
      print("ABORTING LINKING")
      sys.exit(1)


  link_obj = "ccs_app: "
  # keep track of link files in a list
  link_deps = []
  # add common files and main
  link_deps.append(config["main"])
  link_deps = link_deps + config["common"]

  # add files that have options
  for k,v in config["options"].items():
    print("parsing config item:", k, v)
    link_file = method_file_mapping[k][v]
    if (link_file != None):
      link_deps.append(link_file)

  # turn array of filenames to a string with object postfix
  link_rule = link_obj + " ".join([x + ".o" for x in link_deps])

  print(link_rule)
  with open("ccs_app.deps", "w") as f:
      f.write(link_rule)

  sys.exit(0)

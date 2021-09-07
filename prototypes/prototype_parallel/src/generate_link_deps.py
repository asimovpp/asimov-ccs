import sys
import os
import json
from collections import OrderedDict

# define must-appear-together sets
in_set = \
[
  {
    "main": "mpi"
  }
]

# define conflicting configurations
conflicts = \
[
]


def pretty_print(d):
    print(json.dumps(d, indent=2))


def check_for_conflicts(config):
  pass
  #  print("checking config for conflicts")
  #  for c in conflicts:
  #    # check for conflicts that appear as subsets: of config
  #    if (c.items() <= config.items()):
  #      print("config has following illegal combination of settings: ", c)
  #      print("ABORTING LINKING")
  #      sys.exit(1)


def get_link_rule(config):
  link_obj = "ccs_app: "
  # keep track of link files in a list
  link_deps = []
  # add common files and main
  link_deps.append(config["main"])
  link_deps = link_deps + config["common"]
  # add files that have options
  link_deps = link_deps + [v for k,v in config["options"].items()]

  # turn array of filenames to a string with object postfix
  return link_obj + " ".join([x + ".o" for x in link_deps])


def apply_config_mapping(config, config_mapping):
  out = {}
  out["main"] = config["main"]
  out["common"] = list(config["common"])
  out["options"] = {}
  opts = out["options"]
  for k,v in config["options"].items():
    m = config_mapping[k][v]
    print("--", k, v, m) 
    if isinstance(m, dict):
      print("  is dict")
      opts.update(m) 
    elif isinstance(m, str):
      print("  is str")
      opts[k] = m
    else:
      raise Exception(m, "matches no config mapping rule")


  return out


if __name__ == "__main__":
  if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

  # get definition of how to map config options to filenames
  with open("config_mapping.json") as f:
    config_mapping = json.load(f)
  
  with open(sys.argv[1]) as f:
    # want to preserve the order of the config file so that overwriting behaviour is clear
    config = json.load(f, object_pairs_hook=OrderedDict)
  print("---config read: ")
  pretty_print(config)

  mapped_config = apply_config_mapping(config, config_mapping)
  print("---mapped config: ")
  pretty_print(mapped_config)
  check_for_conflicts(mapped_config)
  link_rule = get_link_rule(mapped_config)

  print(link_rule)
  with open("ccs_app.deps", "w") as f:
      f.write(link_rule)

  sys.exit(0)

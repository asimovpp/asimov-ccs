import sys
import os
import json
import logging as log
from collections import OrderedDict
import process_dependencies as pdeps

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
    return json.dumps(d, indent=2)


def check_for_conflicts(config):
  pass
  #  print("checking config for conflicts")
  #  for c in conflicts:
  #    # check for conflicts that appear as subsets: of config
  #    if (c.items() <= config.items()):
  #      print("config has following illegal combination of settings: ", c)
  #      print("ABORTING LINKING")
  #      sys.exit(1)


def get_link_rule(config, deps):
  link_obj = "ccs_app: "
  # keep track of link files in a list
  link_deps = []
  # add main and common files
  link_deps.append(config["main"])
  commons = pdeps.find_commons(deps)
  log.debug("commons: %s", " ".join(commons))
  link_deps = link_deps + commons 
  # add files that have options
  link_deps = link_deps + [v for k,v in config["options"].items()]

  # turn array of filenames to a string with object postfix
  return link_obj + " ".join(["obj/" + os.path.basename(x) + ".o" for x in link_deps])


def apply_config_mapping(config, config_mapping):
  out = {}
  out["main"] = config["main"]
  out["options"] = {}
  opts = out["options"]
  for k,v in config["options"].items():
    m = config_mapping[k][v]
    log.debug("-- %s %s %s",k, v, m) 
    if isinstance(m, dict):
      log.debug("  is dict")
      opts.update(m) 
    elif isinstance(m, str):
      log.debug("  is str")
      opts[k] = m
    else:
      raise Exception(m, "matches no config mapping rule")


  return out


if __name__ == "__main__":
  if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")
  log.basicConfig(stream=sys.stderr, level=log.INFO)

  # get definition of how to map config options to filenames
  # assume the mapping file is in the same directory as this script
  with open(sys.path[0] + "/config_mapping.json") as f:
    config_mapping = json.load(f)
  
  with open(sys.argv[1]) as f:
    # want to preserve the order of the config file so that overwriting behaviour is clear
    config = json.load(f, object_pairs_hook=OrderedDict)
  log.debug("config read:\n%s", pretty_print(config))

  deps = pdeps.parse_dependencies(sys.argv[2])

  mapped_config = apply_config_mapping(config, config_mapping)
  log.debug("mapped config:\n%s", pretty_print(mapped_config))
  check_for_conflicts(mapped_config)
  link_rule = get_link_rule(mapped_config, deps)

  log.info("Configurator produced link rule:\n%s", link_rule)
  with open(sys.argv[3], "w") as f:
      f.write(link_rule)

  sys.exit(0)

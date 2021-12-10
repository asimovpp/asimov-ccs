import sys
import os
import yaml
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
    return yaml.dump(d)


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
  commons = pdeps.find_commons_custom(deps)
  log.debug("commons: %s", " ".join(commons))
  link_deps = link_deps + commons 
  # add files that have options
  link_deps = link_deps + [v for k,v in config["config"].items()]

  # turn array of filenames to a string with object postfix
  return link_obj + " ".join(["obj/" + os.path.basename(x) + ".o" for x in link_deps])


def apply_config_mapping(config, config_mapping):
  out = {"main": "", "config": {}}
  
  if "main" in config:
    out["main"] = config["main"]
  else:
    raise Exception("config has to specify a main")
  
  if "base" in config:
    base = config["base"]
    if base in config_mapping["bases"]:
        out["config"].update(config_mapping["bases"][base]["defaults"])
    else:
        raise Exception("base not found in config mapping", base)
  else:
    raise Exception("config has to specify a base")

  if "options" in config:
    opts = config_mapping["bases"][base]["options"]
    if opts:
      for k,v in config["options"].items():
        if k in opts:
          if v in opts[k]:
            out["config"][k] = v 
          else:
            raise Exception("option choice not found in config mapping", k)
        else:
          raise Exception("option category not found in config mapping", k)
    else:
      raise Exception("base does not allow options", base)
  
  return out


if __name__ == "__main__":
  if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")
  log.basicConfig(stream=sys.stderr, level=log.INFO)

  # get definition of how to map config options to filenames
  # assume the mapping file is in the same directory as this script
  with open(sys.path[0] + "/config_mapping.yaml") as f:
    config_mapping = yaml.load(f, Loader=yaml.FullLoader)
  
  with open(sys.argv[1]) as f:
    # TODO: make sure yaml dictionaries are loaded in the same order they are written
    # want to preserve the order of the config file so that overwriting behaviour is clear
    config = yaml.load(f, Loader=yaml.FullLoader)
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

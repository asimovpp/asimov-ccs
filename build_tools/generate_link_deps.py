import sys
import os
import logging as log
import yaml
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


def get_link_rule(link_obj, config, deps):
  # keep track of link files in a list
  link_deps = []
  # add main and common files
  link_deps.append(config["main"])
  commons = pdeps.find_commons_custom(deps)
  log.debug("commons: %s", " ".join(commons))
  link_deps = link_deps + commons
  # add files that have options
  link_deps = link_deps + [v for k,v in config["config"].items()]
  add_config_extras(link_deps, config)

  # turn array of filenames to a string with object postfix
  return link_obj + ": "+ " ".join(["${OBJ_DIR}/" + os.path.basename(x) + ".o" for x in link_deps])


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

  if "extra" in config:
    out["extra"] = config["extra"]

  return out


def merge_dep_graphs(a, b):
  c = {}
  c.update(a)
  for k,v in b.items():
    if k not in c:
      c[k] = v
    else:
      c[k] = list(set(c[k] + v))
  return c


def generate_minimal_deps(alldeps, main, submods_filename, config):
  # minimise the dependnecy graph with modules needed for main
  deps = pdeps.minimise_deps(alldeps, main)

  # add required submodules and their dependencies to the graph
  # this may add new modules to the graph, so iterate until convergence
  submods = pdeps.parse_submodules(submods_filename)
  psmods = pdeps.find_possible_submods(deps, submods)
  old_psmods = []
  while len(psmods) != len(old_psmods):
    for smod in psmods:
      if smod in config["config"].values():
        smod_deps = pdeps.minimise_deps(alldeps, smod)
        deps = merge_dep_graphs(deps, smod_deps)
    old_psmods = psmods
    psmods = pdeps.find_possible_submods(deps, submods)

  return deps


def add_config_extras(deps, config):
  if "extra" in config:
    extra = config["extra"]
    if isinstance(extra, list):
      for e in extra:
        if isinstance(deps, dict):
          deps[e] = []
        else:
          deps.append(e)
    else:
      if isinstance(deps, dict):
        deps[extra] = []
      else:
        deps.append(extra)


def get_min_link_rule(link_obj, mindeps):
  return link_obj + ": " + " ".join(["${OBJ_DIR}/" + os.path.basename(x) + ".o" for x in mindeps.keys()])


if __name__ == "__main__":
  if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")
  log.basicConfig(stream=sys.stderr, level=log.INFO)

  # get definition of how to map config options to filenames
  # assume the mapping file is in the same directory as this script
  if os.environ['CCS_PROPRIETARY'] == 'yes':
    config_mapping_filename = os.environ['CCS_PROPRIETARY_DIR'] + "/src/build_tools/config_mapping.yaml"
  else:
    config_mapping_filename = sys.path[0] + "/config_mapping.yaml"

  with open(config_mapping_filename) as f:
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

  if len(sys.argv) > 4:
    mindeps = generate_minimal_deps(deps, mapped_config["main"], sys.argv[4], mapped_config)
    add_config_extras(mindeps, mapped_config)
    link_rule = get_min_link_rule("ccs_app", mindeps)
    lib_rule = get_min_link_rule("lib", mindeps)
  else:
    link_rule = get_link_rule("ccs_app", mapped_config, deps)
    link_rule = get_link_rule("lib", mapped_config, deps)

  log.debug("Configurator produced link rule:\n%s", link_rule)
  with open(sys.argv[3], "w") as f:
    f.write(link_rule)
    f.write("\n")
    f.write(lib_rule)

  sys.exit(0)

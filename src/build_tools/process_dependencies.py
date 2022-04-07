import sys, os

def trim(s):
  return s[:-2]


def parse_dependencies(filename):
  data = {}
  with open(filename, "r") as f:
    for line in f:
      obj_files = [trim(x) for x in line.split(" ") if x[-2:] == ".o"]   
      target = os.path.basename(obj_files[0])
      deps = [os.path.basename(x) for x in obj_files[1:]]
      if target in data:
        raise Exception("Error: found duplicate filename in source files:", target)
      data[target] = deps 
  return data


def parse_submodules(filename):
  submods = parse_dependencies(filename)
  for k,v in submods.items():
    if len(v) == 1:
      submods[k] = v[0]
    else:
      raise Exception("Unexpected number of elements encountered when parsing submodules")
  return submods


def find_possible_submods(deps, submods):
  psub = []
  mods = set(deps.keys())
  for smod,mod in submods.items():
    if mod in mods:
      psub.append(smod)

  return set(psub)


def draw_dependencies(deps):
  import matplotlib.pyplot as plt
  import networkx as nx
  g = nx.DiGraph(deps).reverse()
  print(g)
  #nx.write_adjlist(g, "deps_adjlist.txt")
  #nx.draw_shell(g, arrowsize=3, with_labels=True)
  pos = nx.shell_layout(g)
  nx.draw(g, pos, arrowsize=10, node_color="yellow")
  nx.draw_networkx_labels(g, pos, bbox=dict(facecolor='yellow', alpha=0.7), horizontalalignment="left", verticalalignment="top")
  plt.show()


def draw_dependencies_interactive(deps):
  from pyvis.network import Network
  import networkx as nx
  nt = Network('1000px', '1000px', directed=True)
  nt.toggle_physics(False)
  
  g = nx.DiGraph(deps).reverse()
  layout = nx.spiral_layout(g)
  #layout = nx.kamada_kawai_layout(g)
  
  nt.from_nx(g)
  scale = 1000
  for node in nt.nodes:
    node["x"] = layout[node["id"]][0] * scale
    node["y"] = layout[node["id"]][1] * scale
  nt.write_html('deps.html')


def find_commons(deps):
  import networkx as nx
  g = nx.DiGraph(deps).reverse()
  # mains and submodules are leaf nodes, so the internal nodes should be common to all builds
  internal_nodes = [x for x in g.nodes() if g.out_degree(x)!=0]
  return internal_nodes


# return a list of nodes that point to one or more other nodes
# i.e. modules that _are used_ by one or more other modules
def find_commons_custom(deps):
  all_inward = [v for k,v in deps.items()]
  all_inward_flattened = [x for l in all_inward for x in l]
  internal_nodes = set(all_inward_flattened)
  return list(internal_nodes)


def minimise_deps(deps, main):
  min_deps = {}
  next_deps = [main]
  while len(next_deps) > 0:
    dep = next_deps.pop()
    #if dep is already in min_deps, we have encountered a loop
    if not dep in min_deps: 
      min_deps[dep] = deps[dep]
      for d in deps[dep]:
        next_deps.append(d)
    else:
      pass

  return min_deps

if __name__ == "__main__":
  deps = parse_dependencies(sys.argv[1])
  main = sys.argv[2]
  submods = parse_submodules(sys.argv[3])
  min_deps = minimise_deps(deps, main, submods)
  draw_dependencies_interactive(min_deps)

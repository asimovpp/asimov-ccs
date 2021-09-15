import sys

def trim(s):
  return s[:-2]


def parse_dependencies(filename):
  data = {}
  with open(filename, "r") as f:
    for line in f:
      obj_files = [x for x in line.split(" ") if x[-2:] == ".o"]   
      data[trim(obj_files[0])] = [trim(x) for x in obj_files[1:]] 
  return data


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
  g = nx.DiGraph(deps).reverse()
  nt = Network('1000px', '1000px', directed=True)
  nt.from_nx(g)
  nt.toggle_physics(False)
  nt.write_html('deps.html')


def find_commons(deps):
  import networkx as nx
  g = nx.DiGraph(deps).reverse()
  # mains and submodules are leaf nodes, so the internal nodes should be common to all builds
  internal_nodes = [x for x in g.nodes() if g.out_degree(x)!=0]
  return internal_nodes


if __name__ == "__main__":
  deps = parse_dependencies(sys.argv[1])
  print(deps)
  print("common:", find_commons(deps))
  #draw_dependencies(deps)
  draw_dependencies_interactive(deps)

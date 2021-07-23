import sys
import os
import json

conflicts = \
[
  {"stepper": "2loop", "particles": "none"}
]

method_file_mapping = \
{
  "stepper"   : {"2loop": "stepper_2loop",
                 "3loop": "stepper_3loop"},
  "turbulence": {"ke": "turbulence_ke",
                 "kw": "turbulence_kw"},
  "particles" : {"1": "particles_1",
                 "2": "particles_2",
                 "none": None},
  "flux"      : {"1st_order": "flux_1order",
                 "2nd_order": "flux_2order"},
  "solver"    : {"amg": "solve_amg",
                 "cgstab": "solve_cgstab"}
}

if __name__ == "__main__":
  with open(sys.argv[1]) as f:
    config = json.load(f)
  print("config read: ", config)

  print("checking config for conflicts")
  for c in conflicts:
    if (c.items() <= config.items()):
      print("config has following illegal combination of settings: ", c)
      print("ABORTING LINKING")
      exit(1)

  link_command = os.environ["FC"] + " -o ccs_main main.o types.o presolve.o flowsolve.o presolve_basic.o" 
  for k,v in config.items():
    print("parsing config item:", k, v)
    link_file = method_file_mapping[k][v]
    if (link_file != None):
      link_command = link_command + " " + link_file + ".o"
  print("link command: ", link_command)

  os.system(link_command) 

  

import sys
import os
import json

method_file_mapping = \
{
  "stepper"   : {"2loop": "stepper_2loop",
                 "3loop": "stepper_3loop"},
  "turbulence": {"ke": "turbulence_ke",
                 "kw": "turbulence_kw"},
  "particles" : {"1": "particles_1",
                 "2": "particles_2"},
  "flux"      : {"1st_order": "flux_1order",
                 "2nd_order": "flux_2order"},
  "solver"    : {"amg": "solve_amg",
                 "cgstab": "solve_cgstab"}
}

if __name__ == "__main__":
  with open(sys.argv[1]) as f:
    config = json.load(f)
  print("config read: ", config)

  link_command = os.environ["FC"] + " -o ccs_main main.o presolve.o postsolve.o flowsolve.o presolve_basic.o postsolve_basic.o" 
  for k,v in config.items():
    print("parsing config item:", method_file_mapping[k][v])
    link_command = link_command + " " + method_file_mapping[k][v] + ".o"
  print("link command: ", link_command)

  os.system(link_command) 

  

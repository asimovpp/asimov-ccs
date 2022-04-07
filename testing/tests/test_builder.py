import sys
import yaml
import os
import subprocess as sp

config_file = sys.argv[1]
output_stub = sys.argv[2]
test_src_obj = {x: output_stub + os.path.splitext(x)[0] + ".o" for x in sys.argv[3:]}

test_deps = output_stub + ".deps"
smod_deps = output_stub + ".smod.deps"
test_link = output_stub + ".link"
test_exe = output_stub 

print("TEST BUILDER: starting")

print("TEST BUILDER: compiling test objects")
for src,obj in test_src_obj.items():
  sp.run("${FC} ${FFLAGS} -I${OBJ_DIR} ${INC} -c -o " + obj + " " + src,
         shell=True, check=True)

print("TEST BUILDER: evaluating test dependencies")
submods_supported     = "makedepf90 -S " + smod_deps + " ${SRC} " + " ".join(test_src_obj.keys()) + " > " + test_deps + " 2> " + test_deps + ".err"
submods_not_supported = "makedepf90"                 + " ${SRC} " + " ".join(test_src_obj.keys()) + " > " + test_deps + " 2> " + test_deps + ".err"
sp.run("if [ ${MAKEDEPF90_SMODS} = 0 ]; then " + submods_supported + "; else " + submods_not_supported + "; fi",
       shell=True, check=True)

print("TEST BUILDER: generating link rule for test")
submods_supported     = "python3 ${TOOLS}generate_link_deps.py " + config_file + " " + test_deps + " " + test_link + " " + smod_deps
submods_not_supported = "python3 ${TOOLS}generate_link_deps.py " + config_file + " " + test_deps + " " + test_link
sp.run("if [ ${MAKEDEPF90_SMODS} = 0 ]; then " + submods_supported + "; else " + submods_not_supported + "; fi",
       shell=True, check=True)

# getting filename of main
with open(config_file, "r") as f:
  config = yaml.load(f, Loader=yaml.FullLoader)
  main = config["main"]

# stripping off target (and maybe test object file) from link rule
with open(test_link, "r") as f:
  link_rule = f.readline().strip().split(" ")
idx = 0
if len(test_src_obj) == 0:
  idx = 1 #main is found in src/obj, e.g. for case_setup tests
else:
  if any([main in x for x in test_src_obj.keys()]):
    idx = 2 #normal test case
  else:
    idx = 1 #special case for case_setup tests that also have extra files for the test

link_rule = " ".join(link_rule[idx:])

print("TEST BUILDER: linking test executable")
sp.run("${FC} ${FFLAGS} -I${OBJ_DIR} ${INC} " + " ".join(test_src_obj.values()) + " " + link_rule + " -o " + test_exe + " ${LIB}",
       shell=True, check=True)

print("TEST BUILDER: success")

import sys
import yaml
import os
import subprocess as sp

config_file = sys.argv[1]
output_stub = sys.argv[2]

test_obj = output_stub + ".test_obj.o"
test_deps = output_stub + ".deps"
test_link = output_stub + ".link"
test_exe = output_stub 

print("TEST BUILDER: starting")

print("TEST BUILDER: getting path to test source")
with open(config_file, "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    main = os.path.dirname(config_file) + "/" + config["main"] + ".f90"

print("TEST BUILDER: compiling test object")
sp.run("${FC} ${FFLAGS} -I${OBJ_DIR} ${INC} -c -o " + test_obj + " " + main,
       shell=True, check=True)

print("TEST BUILDER: evaluating test dependencies")
sp.run("makedepf90 ${SRC} " + main + " > " + test_deps + " 2> " + test_deps + ".err",
       shell=True, check=True)

print("TEST BUILDER: generating link rule for test")
sp.run("python3 ${TOOLS}generate_link_deps.py " + config_file + " " + test_deps + " " + test_link + " " + os.path.splitext(main)[0] + " submodules.txt",
       shell=True, check=True)

print("TEST BUILDER: stripping off target and test object file from link rule")
with open(test_link, "r") as f:
    link_rule = f.readline().strip()
link_rule = " ".join(link_rule.split(" ")[2:])

print("TEST BUILDER: linking test executable")
sp.run("${FC} ${FFLAGS} -I${OBJ_DIR} ${INC} " + test_obj + " " + link_rule + " -o " + test_exe + " ${LIB}",
       shell=True, check=True)

print("TEST BUILDER: success")

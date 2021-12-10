import sys
import yaml
import os

config_file = sys.argv[1]
test_obj = sys.argv[2]
test_deps = sys.argv[3]
test_link = sys.argv[4]
test_exe = sys.argv[5]

# get path to test source
with open(config_file, "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    main = os.path.dirname(config_file) + "/" + config["main"] + ".f90"

# compile test object
os.system("${FC} ${FFLAGS} -I${OBJ_DIR} ${INC} -c -o " 
           + test_obj + " " + main)

# evaluate test dependencies
os.system("makedepf90 ${SRC} " + main + " > " + test_deps)

# generate link rule for test
os.system("python3 ${TOOLS}generate_link_deps.py " + config_file + " " + test_deps + " " + test_link)

# strip off target and test object file from link rule
with open(test_link, "r") as f:
    link_rule = f.readline().strip()
link_rule = " ".join(link_rule.split(" ")[2:])

# link test executable 
os.system("${FC} ${FFLAGS} -I${OBJ_DIR} ${INC} " 
           + test_obj + " " + link_rule + " -o " + test_exe + " ${LIB}")

"""
This script takes two arguments:
    - inputfile: the yaml config file with CCS configuration
    - extra_to_test: a valid config.yaml "options" value
                     (see build_tools/config_mapping.yaml)
This script will then see if the `extra_to_test` value exists in the
`inputfile` passed and return:
    int: 0 - if `extra_to_test` is not found in `inputfile`
    int: 1 - if `extra_to_test` is found in `inputfile`
"""

import sys
import yaml


if __name__=="__main__":
  inputfile = sys.argv[1]
  extra_to_test = sys.argv[2]
        
  assert isinstance(extra_to_test, str)

  with open(inputfile, "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

  result = 0

  if "options" in config:
    # The first check to avoid ValueErrors
    if extra_to_test in config["options"].values():
      result = 1

  print(result)

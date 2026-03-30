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
    if extra_to_test in config["options"].values():
      result = 1

  print(result)

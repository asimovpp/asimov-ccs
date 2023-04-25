ACTION=$1
TARGET=$2
GOAL_SCORE=8


if [ $ACTION = "lint" ]; then
  for f in $(find $TARGET -name "*.f90"); do
    echo $f
    flint lint -r .flint_config.yaml $f
  done
elif [ $ACTION = "score" ]; then
  score=$(flint score -r .flint_config.yaml $TARGET | grep -Eo '[+-]?[0-9]+[.][0-9]+')
  echo "flint score: "$score"; goal is: "$GOAL_SCORE
  # the comparison sign is the opposite of what one might expect because
  # we need an exit status of 0 to indicate success, but bc returns 0 for False
  exit $(echo "$score < $GOAL_SCORE" | bc -l)
elif [ $ACTION = "score_verbose" ]; then
  flint score -r .flint_config.yaml $TARGET --verbose
elif [ $ACTION = "score_each_file" ]; then
  all_pass=0
  for f in $(find $TARGET -name "*.f90"); do
    echo $f
    score=$(flint score -r .flint_config.yaml $f | grep -Eo '[+-]?[0-9]+[.][0-9]+')
    echo "  flint score: "$score"; goal is: "$GOAL_SCORE
    # returns 0 for False
    if [ $(echo "$score < $GOAL_SCORE" | bc -l) = 1 ]; then
      #check for presence of @dont_fail_linter tag
      #ignore linter score if the tag is found
      if grep -rq "@dont_fail_linter" $f; then
        echo "***FAIL BUT IGNORED***"
      else
        echo "***FAIL***"
        all_pass=1
      fi
      
    fi
  done

  exit $(echo $all_pass) 
elif [ $ACTION = "fprettify" ]; then
  if [ -d $TARGET ]; then
    fprettify -r $TARGET
  else
    fprettify $TARGET
  fi

  # remove spaces around '%'
  find $TARGET -name "*.f90" -exec sed -i 's/ % /%/g' {} \;
  # add spaces around '//'
  find $TARGET -name "*.f90" -exec sed -i 's/\/\// \/\/ /g' {} \;
  if grep -rq "module procedure" $TARGET; then
    echo 'DOUBLE CHECK THAT "module procedure" INDENTATION IS OK IN FILES:'
    grep -rl "module procedure" $TARGET
  fi
else
  echo 'Unknown ACTION provided to linting script'
fi

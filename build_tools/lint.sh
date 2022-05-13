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
fi

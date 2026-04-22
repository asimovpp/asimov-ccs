#!/bin/bash

MPIRUN=$1
#$MPIRUN -n 4 $2 &
$2 &
ID=$!
sleep 2
kill -s TERM $ID
wait $ID
exit_code=$?
echo exit code: $exit_code
exit $exit_code

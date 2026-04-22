#!/bin/bash

$1 -n 4 $2 &
# $2 &
ID=$!
sleep 2
kill -s TERM $ID
wait $ID
exit $?

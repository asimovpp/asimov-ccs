#!/bin/bash

# This test launches a binary that creates a signal trap, then does 
# 2,500,000 loops of busy work. This script waits 2 seconds after the binary
# is launched then sends a sigterm to the test process ID and checks the exit
# code of the test binary. If the binary is allowed to finish the "busy" loop,
# it exits with an error code. If the SEGTERM is caught, it exits with code 0.
# This code is then propagated to the test script, so that `lit` can capture the
# status of the job and tally pass/fail accordingly.

# Run test_sigterm binary
$1 &
# Capture process ID from test binary
ID=$!
sleep 2
# send SIGTERM to process ID running test binary
kill -s TERM $ID
# ensure the test binary has stopped running
wait $ID
# get exit code from test binary and exit script with same code.
exit_code=$?
echo exit code: $exit_code
exit $exit_code

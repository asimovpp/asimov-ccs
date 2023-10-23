#!/bin/bash

# test if average load over last minute is >= 1, i.e. if more than 1 core is doing something
if [ $(echo $(awk '{print $1}' /proc/loadavg) '>= 1.0' | bc -l) ] then
    echo "Load detected on the system:"
    ps -eo pid,user,%cpu,cmd --sort -%cpu | head -n 10
    exit 1
fi

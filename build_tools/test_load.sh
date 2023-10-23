#!/bin/bash

exit $(echo $(awk '{print $1}' /proc/loadavg) '>= 1.0' | bc -l)

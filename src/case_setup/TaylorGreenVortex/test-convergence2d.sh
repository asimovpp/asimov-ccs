#!/bin/bash

export NP=8
export CASE=TaylorGreenVortex2D

export RTOL=1.0e-50
export ATOL=1.0e-16

export CCS=../../../ccs_app

for m in 16 32 64 128 256
do
	mpirun -np ${NP} ${CCS} --ccs_case ${CASE} --ccs_m ${m} -ksp_atol ${ATOL} -ksp_rtol ${RTOL}
	mv tgv2d-err.log tgv2d-err.${m}.log
done

python plot-tgv2derr.py

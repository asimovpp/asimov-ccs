#!/bin/bash

export NP=8
export ATOL=1.0e-16
export RTOL=1.0e-50

for m in 16 32 64 128 256
do
	mpirun -np ${NP} ../../../ccs_app --ccs_case TaylorGreenVortex2D --ccs_m ${m} -ksp_atol ${ATOL} -ksp_rtol ${RTOL}
	mv tgv2d-err.log tgv2d-err.log.${m}
done

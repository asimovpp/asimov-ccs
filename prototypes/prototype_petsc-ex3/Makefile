#
#        FILE: Makefile
# DESCRIPTION: Top-level Makefile for PETSc ex3 prototype of ASiMoV-CCS.
#

MK = make

all: docs ex3

ex3:
	$(MK) -C src

docs: doxy docpdf
docpdf:
	$(MK) -C latex
doxy:
	doxygen doxygen.cfg

clean:
	$(MK) -C src clean
	$(MK) -C latex clean

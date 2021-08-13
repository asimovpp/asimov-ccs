#
#        FILE: Makefile
# DESCRIPTION: Builds/runs tests/documentation for ASiMoV-CCS
#

DOCDIR = docs
TESTDIR = tests

MK = make

all : docs tests

.PHONY : docs
docs :
	$(MK) -C $(DOCDIR)

.PHONY : tests
tests :
	$(MK) -C $(TESTDIR)

clean : clean-docs clean-tests
clean-docs :
	$(MK) -C $(DOCDIR) clean
clean-tests :
	$(MK) -C $(TESTDIR) clean

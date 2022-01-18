#
#        FILE: Makefile
# DESCRIPTION: Builds/runs tests/documentation for ASiMoV-CCS
#

# NOTE: Do not invoke this directly. Call `make tests` from the `src` directory. 

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

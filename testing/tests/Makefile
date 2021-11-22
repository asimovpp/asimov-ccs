#
#        FILE: Makefile
# DESCRIPTION: Builds/runs test cases for ASiMoV-CCS
#

# Test constants
TESTN = 100       # How many random inputs should a test generate?
TESTSCALE = 1.0e9 # How large to scale test inputs? This should cover any conceivable speed (past
                  # light-speed!) and pressures achieved by compressors with PR > 40 and 10^5 Pa
                  # input.


TEST_SUITES =
TEST_SUITES += linear_solvers
TEST_SUITES += mesh

.PHONY: $(TEST_SUITES)
all : $(TEST_SUITES)
	$(foreach suite,$^,$(MK) -C $(suite) all && mpirun $(suite)/$(suite);)

.PHONY: clean
clean : $(TEST_SUITES)
	$(foreach suite,$^,$(MK) -C $(suite) clean;)

MK = make

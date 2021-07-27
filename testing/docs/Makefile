#
#        FILE: Makefile
# DESCRIPTION: Builds documentation for ASiMoV-CCS tests
#

LTX = pdflatex

all : testing.pdf tidy

testing.pdf : testing.tex
	$(LTX) $<
	$(LTX) $<

tidy :
	rm testing.aux testing.log testing.out testing.toc
clean : tidy
	rm testing.pdf 

RUNTESTS = .

all:
	python3 run_lit.py -j 1 -v $(RUNTESTS) --xunit-xml-output="report.xml"
clean:
	rm -rf */Output/ report.xml

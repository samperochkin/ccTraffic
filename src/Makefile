all: ccTraffic.dll
	# comment

ccTraffic.dll: ccTraffic.cpp
	Rscript --vanilla -e "TMB::compile('ccTraffic.cpp','-O0 -g')"
# anothertemplate.so: anothertemplate.cpp
	# Rscript --vanilla -e "TMB::compile('anothertemplate.cpp','-O0 -g')"

clean:
	rm -rf *o

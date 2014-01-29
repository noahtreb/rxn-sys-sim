cpp = g++

cppFlags = -std=c++11 -fopenmp $(dFlags)
cppLibs = -lnetcdf -lnetcdf_c++
oFlags = 
dFlags = -g -Wall -Wextra

srcs = main.cpp
makeObjectFile = $(cpp) $(cppFlags) -c #-o $@

.phony: clean
clean:
	rm -rf FilePrep debug

all: FilePrep

debug: $(srcs)

FilePrep: $(srcs)
	$(cpp) $(cppFlags) -o $@ $(srcs) $(cppLibs)

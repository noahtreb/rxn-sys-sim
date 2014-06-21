# Set C++ compiler.
cpp = g++
#cpp = /home/kaznessi/trebeshn/Software/gcc-4.8.1/bin/g++

# Set the lib and include directories into variables.
libDir = /usr/lib
incDir = /usr/include
#libDir = /home/kaznessi/trebeshn/Software/netcdf-4.3.0/lib
#incDir = /home/kaznessi/trebeshn/Software/netcdf-4.3.0/include

# Set C++ flag.
cppFlags = -I$(incDir) -std=c++11 -fopenmp -pthread $(dFlags)
cppLibs = -L$(libDir) -lnetcdf -lnetcdf_c++
oFlags = -Ofast
dFlags = -g -Wall -Wextra

srcs = main.cpp Distribution.cpp FileInterface.cpp System.cpp PriorityQueue.cpp Reaction.cpp Species.cpp
makeObjectFile = $(cpp) $(cppFlags) -c #-o $@

.phony: clean
clean:
	rm -rf *.o RxnSysSim debug html latex

doc: *.cpp
	doxygen doxygenConfig

all: RxnSysSim

debug: $(srcs)

RxnSysSim: $(srcs)
	$(cpp) $(cppFlags) -o $@ $(srcs) $(cppLibs)

main.o: main.h main.cpp
	$(makeObjectFile) main.cpp

Distribution.o: Distribution.h Distribution.cpp
	$(makeObjectFile) Distribution.cpp 

FileInterface.o: FileInterface.h FileInterface.cpp
	$(makeObjectFile) FileInterface.cpp

System.o: System.h System.cpp
	$(makeObjectFile) System.cpp

PriorityQueue.o: PriorityQueue.h PriorityQueue.cpp
	$(makeObjectFile) PriorityQueue.cpp	

Reaction.o: Reaction.h Reaction.cpp
	$(makeObjectFile) Reaction.cpp

Species.o: Species.h Species.cpp
	$(makeObjectFile) Species.cpp

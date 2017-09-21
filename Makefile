# Makefile for zprime-top-delphes-analysis
# declan.millar@cern.ch

OBJ = get-parameter.o solve-poly.o main.o highest-pt.o match-bjets-to-leps.o neutrino-weighter.o kinematic-reconstructer.o semilepton-reconstructer.o atlas-style.o trim.o progress-bar.o bool-to-string.o analysis.o
BIN = analysis
LIB = lib
SRC = src
OUT = .

ROOTINC = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lTree

HOSTNAME := $(shell hostname)
ifeq ($(HOSTNAME), Lorkhan)
	C = clang++
	BOOSTINC = -I/usr/local/opt/llvm/include -isystem /usr/local/Cellar/boost/1.65.1/include
	BOOSTLIB = "-L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib" -L/usr/local/Cellar/boost/1.65.1/lib
	DELPHESINC = -isystem /Users/declan/Projects/delphes/install/include
	DELPHESLIB = -L /Users/declan/Projects/delphes/install/lib
else ifeq ($(HOSTNAME), cyan03)
	C = g++
	BOOSTINC = -I /local/software/boost/1.61.0/include
	BOOSTLIB = -L /local/software/boost/1.61.0/lib
	DELPHESINC = -I /home/dam1g09/delphes/install/include
	DELPHESLIB = -L /home/dam1g09/delphes/install/lib
else
	C = g++
	BOOSTINC = -I /afs/cern.ch/user/d/demillar/boost_1_64_0
	BOOSTLIB = -L /afs/cern.ch/user/d/demillar/boost_1_64_0/bin.v2/libs
	# BOOSTINC = -I /afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_55
	# BOOSTLIB = -L /afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt/lib
	DELPHESINC = -I /afs/cern.ch/user/d/demillar/delphes/install/include
	DELPHESLIB = -L /afs/cern.ch/user/d/demillar/delphes/install/lib
endif

BOOSTLIBS = $(BOOSTLIB) -lboost_system -lboost_program_options -lboost_filesystem
DELPHESLIBS = $(DELPHESLIB) -lDelphes


CPPFLAGS = $(ROOTINC) $(BOOSTINC) $(DELPHESINC) -Iinclude
LDFLAGS = $(ROOTLIBS) $(BOOSTLIBS) $(DELPHESLIBS)

# Compile all files ending in .cpp in SRC
$(LIB)/%.o: $(SRC)/%.cpp
	$(C) $(CPPFLAGS) -c -o  $@ $<

# Link main file and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(C) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)

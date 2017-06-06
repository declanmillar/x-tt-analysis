# Makefile for zprime-top-delphes-analysis
# declan.millar@cern.ch

OBJ = solve-poly.o main.o two-highest.o match-bjets-to-leps.o neutrino-weighter.o kinematic-reconstructer.o semilepton-reconstructer.o atlas-style.o trim.o progress-bar.o bool-to-string.o analysis.o
BIN = analysis
LIB = lib
SRC = src
OUT = .

ROOTINC = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

HOSTNAME := $(shell hostname)
ifeq ($(HOSTNAME), Sunder)
	BOOSTINC = -isystem /usr/local/Cellar/boost/1.64.0_1/include
	BOOSTLIB = -L /usr/local/Cellar/boost/1.64.0_1/lib
	DELPHESINC = -isystem /Users/declan/Code/delphes/install/include
	DELPHESLIB = -L /Users/declan/Code/delphes/install/lib
else
	# BOOSTINC = -I /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Boost/1.62.0/x86_64-slc6-gcc62-opt/include/boost-1_62
	# BOOSTLIB = -L /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Boost/1.62.0/x86_64-slc6-gcc62-opt/lib
	BOOSTLIB = /afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt/lib
    BOOSTINC = /afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_55
	DELPHESINC = -I /afs/cern.ch/user/d/demillar/delphes/install/include
	DELPHESLIB = -L /afs/cern.ch/user/d/demillar/delphes/install/lib
endif

BOOSTLIBS = $(BOOSTLIB) -lboost_system -lboost_program_options -lboost_filesystem
DELPHESLIBS = $(DELPHESLIB) -lDelphes

C = g++
INC = $(ROOTINC) $(BOOSTINC) $(DELPHESINC) -Iinclude
LIBS = $(ROOTLIBS) $(BOOSTLIBS) $(DELPHESLIBS)

# Compile all files ending in .cpp in SRC
$(LIB)/%.o: $(SRC)/%.cpp
	$(C) $(INC) -c -o  $@ $<

# Link main file and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(C) $(LIBS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)

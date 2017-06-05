# Makefile for zprime-top-analysis
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
else ifeq ($(HOSTNAME), cyan03)
	BOOSTINC = -I /local/software/boost/1.60.0/include
	BOOSTLIB = -L /local/software/boost/1.60.0/lib
else
	BOOSTINC = -I /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Boost/1.62.0/x86_64-slc6-gcc62-opt/include/boost-1_62
	BOOSTLIB = -L /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Boost/1.62.0/x86_64-slc6-gcc62-opt/lib
	DELPHESINC = -I /afs/cern.ch/user/d/demillar/delphes/install/include
	DELPHESLIB = -L /afs/cern.ch/user/d/demillar/delphes/install/lib
endif

BOOSTLIBS = $(BOOSTLIB) -lboost_system -lboost_program_options -lboost_filesystem
DELPHESLIBS = $(DELPHESLIB) -ldelphes

C = c++
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

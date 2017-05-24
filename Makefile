# Makefile for zprime-top-analysis
# declan.millar@cern.ch

OBJ = solve-poly.o main.o neutrino-weighting.o atlas-style.o trim.o progress-bar.o bool-to-string.o analysis.o
BIN = analysis
LIB = lib
SRC = src
OUT = .

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

HOSTNAME := $(shell hostname)
ifeq ($(HOSTNAME), Sunder)
	BOOSTCFLAGS = -isystem /usr/local/Cellar/boost/1.64.0_1/include
	BOOSTLIB = -L /usr/local/Cellar/boost/1.64.0_1/lib
	DELPHESCFLAGS = -isystem /Users/declan/Code/delphes/install/include
	DELPHESLIB = -L /Users/declan/Code/delphes/install/lib
	# DELPHESCFLAGS = -isystem /usr/local/Cellar/madgraph5_amcatnlo/2.5.2/Delphes
	# DELPHESLIB = -L /usr/local/Cellar/madgraph5_amcatnlo/2.5.2/Delphes
else ifeq ($(HOSTNAME), cyan03)
	BOOSTCFLAGS = -I /local/software/boost/1.60.0/include
	BOOSTLIB = -L /local/software/boost/1.60.0/lib
else
	BOOSTCFLAGS = -I /cvmfs/sft.cern.ch/lcg/releases/LCG_87/Boost/1.62.0/x86_64-slc6-gcc49-opt/include
	BOOSTLIB = -L /cvmfs/sft.cern.ch/lcg/releases/LCG_87/Boost/1.62.0/x86_64-slc6-gcc49-opt/lib
	DELPHESCFLAGS = -I /afs/cern.ch/user/d/demillar/delphes/install/include
	DELPHESLIB = -L /afs/cern.ch/user/d/demillar/delphes/install/lib
endif

BOOSTLIBS = $(BOOSTLIB) -lboost_system -lboost_program_options
DELPHESLIBS = $(DELPHESLIB) -ldelphes
# -lExRootAnalysis -lPhysics

C = c++
CFLAGS = $(ROOTCFLAGS) $(BOOSTCFLAGS) $(DELPHESCFLAGS)
LIBS = $(ROOTLIBS) $(BOOSTLIBS) $(DELPHESLIBS)

# Compile all files ending in .cpp in SRC
$(LIB)/%.o: $(SRC)/%.cpp
	$(C) $(CFLAGS) -c -o  $@ $<

# Link main file and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(C) $(LIBS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)

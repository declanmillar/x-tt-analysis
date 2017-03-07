# Makefile for zprime-top-analysis
# declan.millar@cern.ch

OBJ = solve-poly.o main.o atlas-style.o root-tuple.o trim.o progress-bar.o bool-to-string.o analysis.o
BIN = analysis
LIB = lib
SRC = src
OUT = .

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

HOSTNAME := $(shell hostname)
ifeq ($(HOSTNAME), Sunder)
	BOOSTCFLAGS = -isystem /usr/local/Cellar/boost/1.63.0/include
	BOOSTLIB = -L /usr/local/Cellar/boost/1.63.0/lib
	DELPHES = -isystem /afs/cern.ch/user/d/demillar/Delphes-3.4.0
	# -isystem /usr/local/Cellar/madgraph5_amcatnlo/2.5.2/Delphes/external
	DELPHESLIBS = -L /afs/cern.ch/user/d/demillar/Delphes-3.4.0 -L /usr/local/Cellar/madgraph5_amcatnlo/2.5.2/ExRootAnalysis
else ifeq ($(HOSTNAME), cyan03)
	BOOSTCFLAGS = -I /local/software/boost/1.60.0/include
	BOOSTLIB = -L /local/software/boost/1.60.0/lib
else
	BOOSTCFLAGS = -I /afs/cern.ch/sw/lcg/external/Boost/1.60.0/include
	BOOSTLIB = -L /afs/cern.ch/sw/lcg/external/Boost/1.60.0/lib
	DELPHES = -isystem /afs/cern.ch/user/d/demillar/delphes/install/include
	DELPHESLIBS = -L /afs/cern.ch/user/d/demillar/delphes/install/lib
endif

BOOSTLIBS = $(BOOSTLIB) -lboost_system -lboost_program_options

C = c++
CFLAGS = $(ROOTCFLAGS) $(BOOSTCFLAGS) $(DELPHES)
LIBS = $(ROOTLIBS) $(BOOSTLIBS) $(DELPHESLIBS)

# Compile all files ending in .cpp in SRC
$(LIB)/%.o: $(SRC)/%.cpp
	$(C) $(CFLAGS) -c -o  $@ $<

# Link main file and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(C) $(LIBS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)

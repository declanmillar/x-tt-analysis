# Makefile for zprime-top-delphes-analysis
# declan.millar@cern.ch

OBJ = atlas-style.o get-parameter.o solve-quartic.o solve-poly.o main.o highest-pt.o match-bjets-to-leps.o neutrino-weighter.o kinematic-reconstructer.o semilepton-reconstructer.o trim.o progress-bar.o bool-to-string.o analysis.o
BIN = analysis
INC = include
LIB = lib
SRC = src
OUT = .
C = g++

HOSTNAME := $(shell hostname)
ifeq ($(HOSTNAME), cyan03)
	BOOSTINC = -I/local/software/boost/1.61.0/include
	BOOSTLIB = -L/local/software/boost/1.61.0/lib
# elifeq
	# BOOSTINC = ../boost
	# BOOSTLIB = ../boost/bin.v2/libs
endif
BOOSTLIBS = $(BOOSTLIB) -lboost_system -lboost_program_options -lboost_filesystem

ROOTINC = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lTree
DELPHESINC = -isystem ../delphes/install/include
DELPHESLIBS = -L../delphes/install/lib -lDelphes

CPPFLAGS = $(BOOSTINC) $(ROOTINC) $(DELPHESINC) -I$(INC)
LDFLAGS = $(BOOSTLIBS) $(ROOTLIBS) $(DELPHESLIBS)

# Compile all files ending in .cpp in SRC
$(LIB)/%.o: $(SRC)/%.cpp
	$(C) $(CPPFLAGS) -c -o  $@ $<

# Link main file and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(C) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)

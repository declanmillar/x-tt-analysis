# Makefile for zprime-top-analysis
# declan.millar@cern.ch

OBJ = main.o atlas_style.o RootTuple.o analysis.o
BIN = analysis
LIB = Library
SRC = Source
OUT = .

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)

HOSTNAME := $(shell hostname)
ifeq ($(HOSTNAME),Sunder)
	BOOSTFLAGS      = -isystem /usr/local/Cellar/boost/1.60.0_1/include
	BOOSTLIBS       = -L /usr/local/Cellar/boost/1.60.0_1/lib -lboost_system -lboost_program_options
endif
ifeq ($(HOSTNAME),cyan03)
	BOOSTFLAGS      = -I /local/software/boost/1.60.0/include
	BOOSTLIBS       = -L /local/software/boost/1.60.0/lib -lboost_system -lboost_program_options
else
	BOOSTFLAGS      = -I /afs/cern.ch/sw/lcg/external/Boost/1.60.0/include
	BOOSTLIBS       = -L /afs/cern.ch/sw/lcg/external/Boost/1.60.0/lib -lboost_system -lboost_program_options
endif

C               = g++
CFLAGS          = -O -Wall -fPIC -ggdb -std=c++11

L               = g++
LIBS            =

CFLAGS          += $(ROOTCFLAGS) $(BOOSTFLAGS)
LIBS            += $(ROOTLIBS) $(BOOSTLIBS)

# Compile all files ending in .cpp in SRC
$(LIB)/%.o: $(SRC)/%.cpp
	$(C) $(CFLAGS) -c -o  $@ $<

# Link mainfile and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(L) $(LIBS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)

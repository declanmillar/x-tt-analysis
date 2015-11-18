# Makefile for zprime-top-analysis
# declan.millar@cern.ch

OBJ = main.o atlas_style.o RootTuple.o analysis.o
BIN = analysis
LIB = Library
SRC = Source
OUT = .

ROOTCFLAGS      = $(shell root-config --cflags) -isystem /usr/local/Cellar/boost/1.59.0/include
ROOTLIBS        = $(shell root-config --libs) -L /usr/local/Cellar/boost/1.59.0/lib -lboost_system

C               = g++
CFLAGS          = -O -Wall -fPIC -ggdb -std=c++11

L               = g++
LIBS            =

CFLAGS          += $(ROOTCFLAGS)
LIBS            += $(ROOTLIBS) $(BOOSTLIBS)

# Compile all files ending in .cpp in SRC
$(LIB)/%.o: $(SRC)/%.cpp
	$(C) $(CFLAGS) -c -o  $@ $<

# Link mainfile and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(L) $(LIBS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)
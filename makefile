OBJ            =  main.o atlas_style.o overlay.o superpose2d.o RootTuple.o analysis.o

BIN = go
LIB = Library
SRC = Source
OUT = .

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)

C               = g++  
CFLAGS          = -O -Wall -fPIC -ggdb

L               = g++
LIBS            =

CFLAGS          += $(ROOTCFLAGS)   
LIBS            += $(ROOTLIBS)

# Compile all files ending in .cpp in SRC
$(LIB)/%.o: $(SRC)/%.cpp 
	$(C) $(CFLAGS) -c -o  $@ $<

# Link mainfile and all processes
$(OUT)/$(BIN): $(patsubst %, $(LIB)/%, $(OBJ))
	$(L) $(LIBS) -o $@ $^

clean:
	rm -f $(LIB)/*.o $(OUT)/$(BIN)
 

CC = g++

PROG_NAME = simulate

OBJ_PATH = obj
_OBJ = topology.o lattice.o main.o
OBJ = $(patsubst %, $(OBJ_PATH)/%, $(_OBJ))

INCLUDE_PATH = ../include

# include all .hpp headers in ../include folder
DEPS = $(INCLUDE_PATH)/*.hpp

LIBS = -lm
CFLAGS = -Wall -I $(INCLUDE_PATH) -std=c++11 -pg

# make objects
$(OBJ_PATH)/%.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# link objects (make executable)
$(PROG_NAME) : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY : test clean cleandata

test :
	./$(PROG_NAME)

clean :
	rm -f $(PROG_NAME) $(OBJ_PATH)/*.o

cleandata :
	rm ../rvsaData/* ../relaxationData/*

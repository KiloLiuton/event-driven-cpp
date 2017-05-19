CC = g++

PROG_NAME = simulate

OBJ_PATH = obj
_OBJ = main.o topology.o lattice.o
OBJ = $(patsubst %, $(OBJ_PATH)/%, $(_OBJ))

PCG_INCLUDE_PATH = include

DEPS = topology.h lattice.h $(PCG_INCLUDE_PATH)/*.hpp

LIBS = -lm
CFLAGS = -Wall -I $(PCG_INCLUDE_PATH) -std=c++11

$(OBJ_PATH)/%.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(PROG_NAME) : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY : test clean

test :
	./main

clean :
	rm -f $(PROG_NAME) $(OBJ_PATH)/*.o

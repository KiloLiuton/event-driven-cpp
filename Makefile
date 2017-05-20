CC = g++

PROG_NAME = simulate

OBJ_PATH = obj
_OBJ = topology.o lattice.o main.o
OBJ = $(patsubst %, $(OBJ_PATH)/%, $(_OBJ))

PCG_INCLUDE_PATH = include

DEPS = topology.h lattice.h $(PCG_INCLUDE_PATH)/*.hpp

LIBS = -lm
CFLAGS = -Wall -I $(PCG_INCLUDE_PATH) -std=c++11

# make objects
$(OBJ_PATH)/%.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# make executable
$(PROG_NAME) : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY : test clean

test :
	./$(PROG_NAME)

clean :
	rm -f $(PROG_NAME) $(OBJ_PATH)/*.o

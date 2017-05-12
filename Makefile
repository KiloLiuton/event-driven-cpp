CC = g++

PROG_NAME = simulate

OBJ_PATH = obj
_OBJ = main.o mylibimpl.o
OBJ = $(patsubst %, $(OBJ_PATH)/%, $(_OBJ))

PCG_INCLUDE_PATH = /home/kevin/PCGrandom/include

DEPS = mylib.hpp $(PCG_INCLUDE_PATH)/*.hpp

LIBS = -lm
CFLAGS = -Wall -I $(PCG_INCLUDE_PATH)

$(OBJ_PATH)/%.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(PROG_NAME) : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY : test clean

test :
	./main

clean :
	rm -f $(PROG_NAME) $(OBJ_PATH)/*.o

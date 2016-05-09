# Makefile for my_needle.cpp
CC         = g++
CFLAGS     = -Wall -O3
DEBUGFLAGS = -Wall -DDEBUG -ggdb
STATICFLAGS = -Wall -O3 -static
SRC        = my_needle.cpp
OBJ        = $(SRC:.cpp=.o) 
DEBUGOBJ   = $(SRC:.cpp=_debug.o) 
STATICOBJ   = $(SRC:.cpp=_static.o) 
EXE        = my_needle
DEBUGEXE   = my_needle_debug
STATICEXE   = my_needle_static
RM         = /bin/rm -f
CP         = /bin/cp -f

#compile and assemble C++/C source files into object files  -c flag tells the compiler to create only OBJ files
$(EXE): $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ)  -o $(EXE)
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c  $(SRC) 

debug: $(DEBUGOBJ)
	$(CC) $(DEBUGFLAGS) $(DEBUGOBJ)  -o $(DEBUGEXE)
%_debug.o:	%.cpp
	$(CC) $(DEBUGFLAGS) -c $< -o $@

static: $(STATICOBJ)
	$(CC) $(STATICFLAGS) $(STATICOBJ)  -o $(STATICEXE)
%_static.o:	%.cpp
	$(CC) $(STATICFLAGS) -c $< -o $@

install:
	$(CP)  $(EXE)  $(BINPATH)

clean:
	$(RM)  $(OBJ) $(DEBUGOBJ)

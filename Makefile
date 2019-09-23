# Makefile

CCC = /usr/bin/c++
CC  = /usr/bin/cc
CFLAGS = -O3 #-pg
HOME = ./
PWD = ./
INCLUDE = -I${PWD}

SRC_SILL = sillburp.cc burp.cc
OBJ_SILL = ${SRC_SILL:.cc=.o} 

.c.o: 
	$(CC) $(CFLAGS) -c  $< $(INCLUDE)
.cc.o: 
	$(CCC) $(CFLAGS) -c  $< $(INCLUDE)

all: sillburp

sillburp: $(OBJ_SILL) sillburp.h burp.h
	$(CCC) $(CFLAGS) -o $@ $(OBJ_SILL) -lm

clean:
	/bin/rm $(OBJ_SILL) sillburp *.time *.xy *.dat

PRGS = MACML 
CC = g++
CFLAGS = -O3

all:	$(PRGS) 

MACML: MACML.cpp HC.cpp base.cpp
	$(CC) $(CFLAGS) -std=c++11 -o  $@ MACML.cpp HC.cpp base.cpp  -w



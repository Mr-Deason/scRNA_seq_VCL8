CC = gcc
CPP = g++
LD = g++ -L /usr/local/lib

LIBS = -lz -lboost_filesystem -lboost_system
OCFLAGS = -std=c++11

get_barcode.o : get_barcode.cpp
	$(LD) $(OCFLAGS) $(LIBS)
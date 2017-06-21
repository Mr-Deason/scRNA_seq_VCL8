CC = gcc
CPP = g++
LD = g++ -L /usr/local/lib

LIBS = -lz -lboost_filesystem -lboost_system
OCFLAGS = -std=c++11

get_barcode : get_barcode.o
	$(LD) -o get_barcode $(OCFLAGS) get_barcode.o $(LIBS)
get_barcode.o : get_barcode.cpp
	$(CPP) $(OCFLAGS) -c get_barcode.cpp
CC = gcc
CPP = g++
LD = g++ -L /usr/local/lib

LIBS = -lz -lboost_filesystem -lboost_system
OCFLAGS = -std=c++11

get_barcode : get_barcode.o
	$(LD) -o get_barcode $(OCFLAGS) get_barcode.o $(LIBS)
get_barcode.o : get_barcode.cpp
	$(CPP) $(OCFLAGS) -c get_barcode.cpp
error_correct_and_split : error_correct_and_split.o
	$(LD) -o error_correct_and_split $(OCFLAGS) error_correct_and_split.o $(LIBS)
error_correct_and_split.o : error_correct_and_split.cpp
	$(CPP) $(OCFLAGS) -c error_correct_and_split.cpp
compute_TCCs : compute_TCCs.o
	$(LD) -o compute_TCCs $(OCFLAGS) compute_TCCs.o $(LIBS)
compute_TCCs.o : compute_TCCs.cpp
	$(CPP) $(OCFLAGS) -c compute_TCCs.cpp
prep_TCC_matrix : prep_TCC_matrix.o
	$(LD) -o prep_TCC_matrix $(OCFLAGS) prep_TCC_matrix.o $(LIBS)
prep_TCC_matrix.o : prep_TCC_matrix.cpp
	$(CPP) $(OCFLAGS) -c prep_TCC_matrix.cpp

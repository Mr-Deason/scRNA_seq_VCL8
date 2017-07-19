CC = gcc
CPP = g++
LD = g++ -L /usr/local/lib

LIBS = -lz -lboost_filesystem -lboost_system
CFLAG = -std=c++11 


default: all
all: get_barcode error_correct_and_split compute_TCCs prep_TCC_matrix gen_sub

get_barcode : get_barcode.o
	$(LD) -o get_barcode $(CFLAG) get_barcode.o $(LIBS)
get_barcode.o : get_barcode.cpp
	$(CPP) $(CFLAG) -c get_barcode.cpp
error_correct_and_split : error_correct_and_split.o
	$(LD) -o error_correct_and_split $(CFLAG) error_correct_and_split.o $(LIBS)
error_correct_and_split.o : error_correct_and_split.cpp
	$(CPP) $(CFLAG) -c error_correct_and_split.cpp
compute_TCCs : compute_TCCs.o
	$(LD) -o compute_TCCs $(CFLAG) compute_TCCs.o $(LIBS)
compute_TCCs.o : compute_TCCs.cpp
	$(CPP) $(CFLAG) -c compute_TCCs.cpp
prep_TCC_matrix : prep_TCC_matrix.o
	$(LD) -o prep_TCC_matrix $(CFLAG) prep_TCC_matrix.o $(LIBS)
prep_TCC_matrix.o : prep_TCC_matrix.cpp
	$(CPP) $(CFLAG) -c prep_TCC_matrix.cpp
gen_sub : gen_sub.o
	$(LD) -o gen_sub $(CFLAG) gen_sub.o $(LIBS)
gen_sub.o : gen_sub.cpp
	$(CPP) $(CFLAG) -c gen_sub.cpp

clean:
	rm get_barcode error_correct_and_split compute_TCCs prep_TCC_matrix gen_sub
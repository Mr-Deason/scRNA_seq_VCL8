CC = gcc
CPP = g++
LD = g++ -L /usr/local/lib

LIBS = -lz -lboost_filesystem -lboost_system -lopenblas -lpthread
CFLAG = -fopenmp -std=c++11
# -O3 -mavx -msse -msse2 -mavx2
INCLU = -I /home/kuangdi/downloads/eigen-3.3.4/


default: all
all: get_barcode \
	error_correct_and_split \
	compute_TCCs \
	prep_TCC_matrix \
	gen_sub \
	gen_shuffle_sub \
	new_error_correct_and_split \
	prep_TCC_matrix_cblas \
	prep_TCC_matrix_eigen \
	prep_TCC_matrix_simd

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
gen_shuffle_sub : gen_shuffle_sub.o
	$(LD) -o gen_shuffle_sub $(CFLAG) gen_shuffle_sub.o $(LIBS)
gen_shuffle_sub.o : gen_shuffle_sub.cpp
	$(CPP) $(CFLAG) -c gen_shuffle_sub.cpp

new_error_correct_and_split : new_error_correct_and_split.o
	$(LD) -o new_error_correct_and_split $(CFLAG) new_error_correct_and_split.o $(LIBS)
new_error_correct_and_split.o : new_error_correct_and_split.cpp
	$(CPP) $(CFLAG) -c new_error_correct_and_split.cpp

prep_TCC_matrix_cblas : prep_TCC_matrix_cblas.o
	$(LD) -o prep_TCC_matrix_cblas $(CFLAG) prep_TCC_matrix_cblas.o $(LIBS)
prep_TCC_matrix_cblas.o : prep_TCC_matrix_cblas.cpp
	$(CPP) $(CFLAG) -c prep_TCC_matrix_cblas.cpp

prep_TCC_matrix_eigen : prep_TCC_matrix_eigen.o
	$(LD) -o prep_TCC_matrix_eigen $(CFLAG) prep_TCC_matrix_eigen.o $(LIBS)
prep_TCC_matrix_eigen.o : prep_TCC_matrix_eigen.cpp
	$(CPP) $(INCLU) $(CFLAG) -c prep_TCC_matrix_eigen.cpp

prep_TCC_matrix_simd : prep_TCC_matrix_simd.o
	$(LD) -o prep_TCC_matrix_simd $(CFLAG) prep_TCC_matrix_simd.o $(LIBS)
prep_TCC_matrix_simd.o : prep_TCC_matrix_simd.cpp
	$(CPP) $(INCLU) $(CFLAG) -c prep_TCC_matrix_simd.cpp

clean:
	rm get_barcode \
	error_correct_and_split \
	compute_TCCs \
	prep_TCC_matrix \
	gen_sub \
	gen_shuffle_sub \
	new_error_correct_and_split \
	prep_TCC_matrix_cblas \
	prep_TCC_matrix_eigen \
	prep_TCC_matrix_simd \
	*.o
# Makefile for SVM-struct, 03.10.06

#Use the following to compile under unix or cygwin
CC = g++
LD = g++

# Added by Al
INCLUDES = -I../../../../core -I../../../../lib/libsvm-2.89/ `pkg-config --cflags opencv`

#Call 'make' using the following line to make CYGWIN produce stand-alone Windows executables
#		make 'SFLAGS=-mno-cygwin'

CFLAGS =    $(SFLAGS) -O3 -fopenmp -fomit-frame-pointer -ffast-math -Wall $(INCLUDES)
LDFLAGS =   $(SFLAGS) -O3 -lm -Wall -fopenmp
#CFLAGS =   $(SFLAGS) -g -pg -Wall $(INCLUDES)
#LDFLAGS =  $(SFLAGS) -g -pg -Wall

all: svm_struct_noexe

svm_struct_noexe: svm_struct_learn.o svm_struct_classify.o svm_struct_common.o svm_struct_main.o svm_struct_main_loss_update.o 

.PHONY: clean
clean:
	rm -f *.o *.tcov *.d core gmon.out *.stackdump


#----------------------#
#----  STRUCT SVM  ----#
#----------------------#

svm_struct_common.o: svm_struct_common.c svm_struct_common.h ../svm_struct_api_types.h
	$(CC) -c $(CFLAGS) svm_struct_common.c -o svm_struct_common.o

svm_struct_learn.o: svm_struct_learn.c ../svm_light/svm_learn.h svm_struct_common.h ../svm_struct_api.h ../svm_struct_api_types.h
	$(CC) -c $(CFLAGS) svm_struct_learn.c -o svm_struct_learn.o

svm_struct_main.o: svm_struct_main.c ../svm_light/svm_common.h ../svm_light/svm_learn.h svm_struct_learn.h svm_struct_common.h ../svm_struct_api.h ../svm_struct_api_types.h
	$(CC) -c $(CFLAGS) svm_struct_main.c -o svm_struct_main.o

svm_struct_classify.o: svm_struct_classify.c svm_struct_common.h ../svm_struct_api_types.h ../svm_struct_api.h ../svm_light/svm_common.h 
	$(CC) -c $(CFLAGS) svm_struct_classify.c -o svm_struct_classify.o

svm_struct_main_loss_update.o: svm_struct_main_loss_update.c ../svm_light/svm_common.h ../svm_light/svm_learn.h svm_struct_learn.h svm_struct_common.h ../svm_struct_api.h ../svm_struct_api_types.h
	$(CC) -c $(CFLAGS) svm_struct_main_loss_update.c -o svm_struct_main_loss_update.o

iterative_svm_struct_main.o: iterative_svm_struct_main.c ../svm_light/svm_common.h ../svm_light/svm_learn.h svm_struct_learn.h svm_struct_common.h ../svm_struct_api.h ../svm_struct_api_types.h
	$(CC) -c $(CFLAGS) iterative_svm_struct_main.c -o iterative_svm_struct_main.o


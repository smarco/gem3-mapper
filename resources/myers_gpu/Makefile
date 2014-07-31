
# Shell interpreter #
#####################
SHELL := /bin/bash

# Compilers  #
##############
NVCC=nvcc
CC=gcc
ICC=icc

# CPU flags  #
##############
CFLAGS=-O3 -m64
LFLAGS=

# CUDA flags #
##############
CUDA_OPT_FLAGS=-O3 -m64
CUDA_LIBRARY_FLAGS=-L/usr/local/cuda/lib64 -I/usr/local/cuda/include -lcuda -lcudart

# DEBUG FLAGS #
###############
GCC_GDB_FLAGS=-g -O0 -m64
NVCC_GDB_FLAGS=-g -G -O0 -m64

ifeq ($(arch),tesla)
	NVCFLAGS=-gencode arch=compute_10,code=sm_10 --ptxas-options=-v
	CUDA_NUM_THREADS=-DCUDA_NUM_THREADS=$(th_block)
	CUDA_NUM_THREADS_PER_SM=-DMAX_THREADS_PER_SM=$(th_sm)
endif

ifeq ($(arch),fermi)
	NVCFLAGS=-gencode arch=compute_20,code=sm_20 --ptxas-options=-v
	CUDA_DEVICE=-DDEVICE=1
	CUDA_NUM_THREADS=-DCUDA_NUM_THREADS=$(th_block)
	CUDA_NUM_THREADS_PER_SM=-DMAX_THREADS_PER_SM=$(th_sm)
endif

ifeq ($(arch),kepler)
	NVCFLAGS=-gencode arch=compute_30,code=sm_30 --ptxas-options=-v
	CUDA_DEVICE=-DDEVICE=0
	CUDA_NUM_THREADS=-DCUDA_NUM_THREADS=$(th_block)
	CUDA_NUM_THREADS_PER_SM=-DMAX_THREADS_PER_SM=$(th_sm)
endif

ifeq ($(arch),titan)
	NVCFLAGS=-gencode arch=compute_35,code=sm_35 --ptxas-options=-v
	CUDA_DEVICE=-DDEVICE=0
	CUDA_NUM_THREADS=-DCUDA_NUM_THREADS=$(th_block)
	CUDA_NUM_THREADS_PER_SM=-DMAX_THREADS_PER_SM=$(th_sm)
endif

# Conf save results #
#####################
MEMORY_WRITES=-DMEMORY_WRITES=1

ifeq ($(funnel), on)
	OPT_FUNNEL=-DFUNNEL
endif

ifeq ($(shuffle), on)
	OPT_SHUFFLE=-DSHUFFLE
endif

ifeq ($(ballot), on)
	OPT_BALLOT=-DBALLOT
endif

all: release

######################
#  BUILD PARSERS
######################

gcandidates-GEM-1st:
	$(CC) $(CFLAGS) genCandidate-reorder.c -o gcandidates-GEM-1st

##################################################
#  MYERS CPU VERSIONS		  					 #
##################################################

myersCPU-padding-32bits:
	$(CC) $(CFLAGS) -o myersCPU-padding-32bits myers-4bits-padding-32bits.c $(LFLAGS) -lrt -fopenmp

myersCPU-GEM-1st:
	$(CC) $(CFLAGS) -o myersCPU-GEM-1st myers-reorderBuffer.c $(LFLAGS) -lrt -fopenmp

myersCPU-GEM-1st-Padding:
	$(CC) $(CFLAGS) -o myersCPU-GEM-1st-Padding myers-reorderBuffer-padding.c $(LFLAGS) -lrt -fopenmp

##################################################
#  1st STEP INTEGRATION		(reorder buffer)	 #
##################################################

release: 
	nvcc -O3 -m64 -gencode arch=compute_20,code=sm_20 --ptxas-options=-v -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	nvcc -O3 -m64 -gencode arch=compute_30,code=sm_30 --ptxas-options=-v -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	nvcc -O3 -m64 -gencode arch=compute_35,code=sm_35 --ptxas-options=-v -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	nvcc -O3 -m64 -gencode arch=compute_50,code=sm_50 --ptxas-options=-v -c coreMyersGPU-CC50.cu -o ../build/coreMyersGPU50.o
	gcc -O3 -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o

debug: 
	nvcc -G -O0 -m64 -gencode arch=compute_20,code=sm_20 --ptxas-options=-v -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	nvcc -G -O0 -m64 -gencode arch=compute_30,code=sm_30 --ptxas-options=-v -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	nvcc -G -O0 -m64 -gencode arch=compute_35,code=sm_35 --ptxas-options=-v -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	nvcc -G -O0 -m64 -gencode arch=compute_50,code=sm_50 --ptxas-options=-v -c coreMyersGPU-CC50.cu -o ../build/coreMyersGPU50.o
	gcc  -g -O0 -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o

clean:
  	
######################
#  INSTALL
######################

install:
	find ./ -perm /a+x -type f -exec mv {} ../bin \;

# Shell interpreter #
#####################
SHELL := /bin/bash

# Compilers  #
##############
NVCC=nvcc
CC=gcc

# CPU flags  #
##############
CFLAGS=-O3 -m64
LFLAGS=

# CUDA flags #
##############
CUDA_FLAGS=-O3 -m64 -gencode arch=compute_35,code=sm_35
CUDA_LIBRARY_FLAGS=-L/usr/local/cuda/lib64 -I/usr/local/cuda/include -lcuda -lcudart

# DEBUG FLAGS #
###############
GCC_GDB_FLAGS=-g -O0 -m64
NVCC_GDB_FLAGS=-g -G -O0 -m64 -gencode arch=compute_52,code=sm_52

all:
	gcc -O3 -c src/gpu_sample.c -o build/gpu_sample.o -I /usr/local/cuda/include/
	gcc -O3 -c src/gpu_fmi.c -o build/gpu_fmi.o -I /usr/local/cuda/include/
	gcc -O3 -c src/gpu_devices.c -o build/gpu_devices.o -I /usr/local/cuda/include/
	gcc -O3 -c src/gpu_buffers.c -o build/gpu_buffers.o -I /usr/local/cuda/include/
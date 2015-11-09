SHELL := /bin/bash

NVCC=nvcc
CC=gcc
CFLAGS=-O3 -m64
CUDA_FLAGS=-O3 -m64 -gencode arch=compute_35,code=sm_35
CUDA_INCLUDES=-I/usr/local/cuda/include
CUDA_LIBS=-I/usr/local/cuda/lib64 -lcuda -cudalrt

GCC_GDB_FLAGS=-g -O0 -m64
NVCC_GDB_FLAGS=-g -G -O0 -m64 -gencode arch=compute_35,code=sm_35

#gcc -c -combine a.c b.c -o c.o

all:
	$(CC) $(CFLAGS) -c src/gpu_sample.c -o build/gpu_sample.o $(CUDA_INCLUDES)
	$(CC) $(CFLAGS) -c src/gpu_fmi.c -o build/gpu_fmi.o $(CUDA_INCLUDES)
	$(CC) $(CFLAGS) -c src/gpu_bpm.c -o build/gpu_bpm.o $(CUDA_INCLUDES)
	$(CC) $(CFLAGS) -c src/gpu_devices.c -o build/gpu_devices.o $(CUDA_INCLUDES)
	$(CC) $(CFLAGS) -c src/gpu_buffers.c -o build/gpu_buffers.o $(CUDA_INCLUDES)
	$(CC) $(CFLAGS) -c src/gpu_index.c -o build/gpu_index.o $(CUDA_INCLUDES)
	$(CC) $(CFLAGS) -c src/gpu_reference.c -o build/gpu_reference.o $(CUDA_INCLUDES)
	$(NVCC) $(CUDA_FLAGS) -c src/gpu_fmi_decode.cu -o build/gpu_fmi_decode.o $(CUDA_INCLUDES)
	$(NVCC) $(CUDA_FLAGS) -c src/gpu_fmi_search.cu -o build/gpu_fmi_search.o $(CUDA_INCLUDES)
	$(NVCC) $(CUDA_FLAGS) -c src/gpu_bpm_filter.cu -o build/gpu_bpm_filter.o $(CUDA_INCLUDES)
	ld -r build/gpu_sample.o build/gpu_fmi.o build/gpu_bpm.o build/gpu_devices.o build/gpu_buffers.o build/gpu_index.o build/gpu_reference.o build/gpu_fmi_decode.o build/gpu_fmi_search.o build/gpu_bpm_filter.o -o build/gem_gpu.o 
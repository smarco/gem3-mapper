SHELL := /bin/bash

NVCC=nvcc
CC=gcc
CFLAGS=-O3 -m64
CUDA_FLAGS=-O3 -m64 -gencode arch=compute_35,code=sm_35
CUDA_INCLUDES=-I/usr/local/cuda/include
CUDA_LIBS=-L/usr/local/cuda/lib64 -lcudart

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
	make bpm_benchmark
	make search_benchmark
	make decode_benchmark


bpm_benchmark:
	$(CC) $(CFLAGS) tools/gpu_benchmark_bpm.c build/gem_gpu.o -o bin/gpu_bpm_bench $(CUDA_INCLUDES) $(CUDA_LIBS) -fopenmp -lrt 
search_benchmark:
	$(CC) $(CFLAGS) tools/gpu_benchmark_search.c build/gem_gpu.o -o bin/gpu_search_bench $(CUDA_INCLUDES) $(CUDA_LIBS) -fopenmp -lrt
decode_benchmark:
	$(CC) $(CFLAGS) tools/gpu_benchmark_decode.c build/gem_gpu.o -o bin/gpu_decode_bench $(CUDA_INCLUDES) $(CUDA_LIBS) -fopenmp -lrt
	

debug:
	$(CC) $(GCC_GDB_FLAGS) -c src/gpu_sample.c -o build/gpu_sample.o $(CUDA_INCLUDES)
	$(CC) $(GCC_GDB_FLAGS) -c src/gpu_fmi.c -o build/gpu_fmi.o $(CUDA_INCLUDES)
	$(CC) $(GCC_GDB_FLAGS) -c src/gpu_bpm.c -o build/gpu_bpm.o $(CUDA_INCLUDES)
	$(CC) $(GCC_GDB_FLAGS) -c src/gpu_devices.c -o build/gpu_devices.o $(CUDA_INCLUDES)
	$(CC) $(GCC_GDB_FLAGS) -c src/gpu_buffers.c -o build/gpu_buffers.o $(CUDA_INCLUDES)
	$(CC) $(GCC_GDB_FLAGS) -c src/gpu_index.c -o build/gpu_index.o $(CUDA_INCLUDES)
	$(CC) $(GCC_GDB_FLAGS) -c src/gpu_reference.c -o build/gpu_reference.o $(CUDA_INCLUDES)
	$(NVCC) $(NVCC_GDB_FLAGS) -c src/gpu_fmi_decode.cu -o build/gpu_fmi_decode.o $(CUDA_INCLUDES)
	$(NVCC) $(NVCC_GDB_FLAGS) -c src/gpu_fmi_search.cu -o build/gpu_fmi_search.o $(CUDA_INCLUDES)
	$(NVCC) $(NVCC_GDB_FLAGS) -c src/gpu_bpm_filter.cu -o build/gpu_bpm_filter.o $(CUDA_INCLUDES)
	ld -r build/gpu_sample.o build/gpu_fmi.o build/gpu_bpm.o build/gpu_devices.o build/gpu_buffers.o build/gpu_index.o build/gpu_reference.o build/gpu_fmi_decode.o build/gpu_fmi_search.o build/gpu_bpm_filter.o -o build/gem_gpu.o 
	make bpm_benchmark_debug
	make search_benchmark_debug
	make decode_benchmark_debug


bpm_benchmark_debug:
	$(CC) $(GCC_GDB_FLAGS) tools/gpu_benchmark_bpm.c build/gem_gpu.o -o bin/gpu_bpm_bench $(CUDA_INCLUDES) $(CUDA_LIBS) -fopenmp -lrt -lstdc++
search_benchmark_debug:
	$(CC) $(GCC_GDB_FLAGS) tools/gpu_benchmark_search.c build/gem_gpu.o -o bin/gpu_search_bench $(CUDA_INCLUDES) $(CUDA_LIBS) -fopenmp -lrt -lstdc++
decode_benchmark_debug:
	$(CC) $(GCC_GDB_FLAGS) tools/gpu_benchmark_decode.c build/gem_gpu.o -o bin/gpu_decode_bench $(CUDA_INCLUDES) $(CUDA_LIBS) -fopenmp -lrt -lstdc++


clean:
	rm -f bin/*
	rm -f build/*.o

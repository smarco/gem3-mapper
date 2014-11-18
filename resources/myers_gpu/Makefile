#============================================================================
# PROJECT: Bit-Parallel Myers GPU
# FILE: Makefile
# DATE: 4/9/2014
# DESCRIPTION: makefile
#============================================================================

include ../../Makefile.mk

SHELL := /bin/bash
NVCC_VERSION = $(shell $(NVCC) --version | grep release | sed 's/.*release //' |  sed 's/,.*//' | sed 's/\.//')
CUDA_LIBRARY_FLAGS = $(INC_CUDA) $(LIBS_CUDA)
FOLDER_BUILD=../build
NVCC_MINIMAL_VERSION = $(shell [[ ($(NVCC_VERSION) -ge "0")  && ($(NVCC_VERSION) -le  "29") ]] && echo "10"; \
                               [[ ($(NVCC_VERSION) -ge "30") && ($(NVCC_VERSION) -le  "49") ]] && echo "30"; \
                               [[ ($(NVCC_VERSION) -ge "50") && ($(NVCC_VERSION) -le  "59") ]] && echo "50"; \
                               [[ ($(NVCC_VERSION) -ge "60") ]] && echo "60";)

$(shell test -d $(FOLDER_BUILD) || mkdir $(FOLDER_BUILD))

ifeq ($(NVCC_MINIMAL_VERSION), 10)
	$(error CUDA SDK version $(NVCC_VERSION) NOT SUPPORTED - (At least CUDA SDK 3.0 is required))
endif

CUDA_MODULES_30 =  coreMyers-GPU20
CUDA_MODULES_50 = $(CUDA_MODULES_30) coreMyers-GPU30 coreMyers-GPU35
CUDA_MODULES_60 = $(CUDA_MODULES_50) coreMyers-GPU50

CUDA_MODULES = $(CUDA_MODULES_$(NVCC_MINIMAL_VERSION))
CUDA_SRCS=$(addsuffix .cu, $(CUDA_MODULES))
CUDA_OBJS=$(addprefix $(FOLDER_BUILD)/, $(CUDA_SRCS:.cu=.o))

MODULES = myersGPU
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))
	

all: release

release: NVCC_COMPILE_FLAGS=-O3 -m64 
release: GCC_COMPILE_FLAGS=-O3 -m64
release: $(CUDA_OBJS) $(OBJS)

profile: NVCC_COMPILE_FLAGS=-O3 -m64 --ptxas-options=-v -lineinfo
profile: GCC_COMPILE_FLAGS=-O3 -m64
profile: $(CUDA_OBJS) $(OBJS)

debug: NVCC_COMPILE_FLAGS=-G -O0 -m64 --ptxas-options=-v
debug: GCC_COMPILE_FLAGS=-g -O0 -m64
debug: $(CUDA_OBJS) $(OBJS)



test: NVCC_COMPILE_FLAGS=-O3 -m64
test: GCC_COMPILE_FLAGS=-O3 -m64
test: $(CUDA_OBJS) $(OBJS)
	$(CC) $(GCC_COMPILE_FLAGS) test.c $(CUDA_OBJS) $(OBJS) -o $(FOLDER_BUILD)/test -fopenmp $(CUDA_LIBRARY_FLAGS) -lrt
	$(CC) $(GCC_COMPILE_FLAGS) genCandidate.c -o $(FOLDER_BUILD)/gcandidates
	$(CC) $(GCC_COMPILE_FLAGS) genReference.c -o $(FOLDER_BUILD)/greference


test-profile: NVCC_COMPILE_FLAGS=-O3 -m64 --ptxas-options=-v -lineinfo
test-profile: GCC_COMPILE_FLAGS=-O3 -m64
test-profile: $(CUDA_OBJS) $(OBJS)
	$(CC) $(GCC_COMPILE_FLAGS) test.c $(CUDA_OBJS) $(OBJS) -o $(FOLDER_BUILD)/test-profile -fopenmp $(CUDA_LIBRARY_FLAGS) -lrt
	$(CC) $(GCC_COMPILE_FLAGS) genCandidate.c -o $(FOLDER_BUILD)/gcandidates
	$(CC) $(GCC_COMPILE_FLAGS) genReference.c -o $(FOLDER_BUILD)/greference


test-debug: NVCC_COMPILE_FLAGS=-O0 -G -m64 --ptxas-options=-v -lineinfo
test-debug: GCC_COMPILE_FLAGS=-O0 -g -m64
test-debug: CUDA_LIBRARY_FLAGS+=-lstdc++
test-debug: $(CUDA_OBJS) $(OBJS)
	$(CC) $(GCC_COMPILE_FLAGS) test.c $(CUDA_OBJS) $(OBJS) -o $(FOLDER_BUILD)/test-debug -fopenmp $(CUDA_LIBRARY_FLAGS) -lrt
	$(CC) $(GCC_COMPILE_FLAGS) genCandidate.c -o $(FOLDER_BUILD)/gcandidates
	$(CC) $(GCC_COMPILE_FLAGS) genReference.c -o $(FOLDER_BUILD)/greference

    
$(FOLDER_BUILD)/coreMyers-GPU20.o: coreMyersGPU-CC20.cu
	$(NVCC) $(NVCC_COMPILE_FLAGS) -gencode arch=compute_20,code=\"sm_20,compute_20\" -c $< -o $@

$(FOLDER_BUILD)/coreMyers-GPU30.o: coreMyersGPU-CC30.cu
	$(NVCC) $(NVCC_COMPILE_FLAGS) -gencode arch=compute_30,code=\"sm_30,compute_30\" -c $< -o $@

$(FOLDER_BUILD)/coreMyers-GPU35.o: coreMyersGPU-CC35.cu
	$(NVCC) $(NVCC_COMPILE_FLAGS) -gencode arch=compute_35,code=\"sm_35,compute_35\" -c $< -o $@

$(FOLDER_BUILD)/coreMyers-GPU50.o: coreMyersGPU-CC50.cu
	$(NVCC) $(NVCC_COMPILE_FLAGS) -gencode arch=compute_50,code=\"sm_50,compute_50\" -c $< -o $@
	
$(OBJS): $(SRCS)
	gcc $(GCC_COMPILE_FLAGS) -c $< -o $@ $(CUDA_LIBRARY_FLAGS)



clean:
	rm -f ../build/*.o

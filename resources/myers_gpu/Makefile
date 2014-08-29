SHELL := /bin/bash

NVCC_VERSION = $(shell nvcc --version | grep release | sed 's/.*release //' |  sed 's/,.*//' | sed 's/\.//')
CUDA_LIBRARY_FLAGS = -L/usr/local/cuda/lib64 -I/usr/local/cuda/include -lcuda -lcudart

NVCC_MINIMAL_VERSION = $(shell [[ ($(NVCC_VERSION) -ge "0")  && ($(NVCC_VERSION) -le  "29") ]] && echo "10"; \
                               [[ ($(NVCC_VERSION) -ge "30") && ($(NVCC_VERSION) -le  "49") ]] && echo "30"; \
                               [[ ($(NVCC_VERSION) -ge "50") && ($(NVCC_VERSION) -le  "59") ]] && echo "50"; \
                               [[ ($(NVCC_VERSION) -ge "60") ]] && echo "60";)

all: NVCC_COMPILE_FLAGS=-O3 -m64 
all: GCC_COMPILE_FLAGS=-O3 -m64
all: release_$(NVCC_MINIMAL_VERSION)

######################
#  BUILD PARSERS
######################

gcandidates:
	gcc -O3 -m64 genCandidate.c -o ../build/gcandidates

gcandidates-debug:
	gcc -O0 -g genCandidate.c -o ../build/gcandidates-debug

greference:
	gcc -O3 -m64 genReference.c -o ../build/greference

greference-debug:
	gcc -O0 -g genReference.c -o ../build/greference-debug

##############################
#  MYERS RELEASE VERSIONS	 #
##############################

#TODO: check compiler version
#TODO: fix the objects requeriments (in the parent makefile)
#INCLUDED FATBINARY COMPILATION TO SUPPORT NEXT ARCHITECTURES

FOLDER_BUILD=../build
MODULES=coreMyersGPU20 coreMyersGPU30 coreMyersGPU35 coreMyersGPU50 myersGPU

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

release_10:
	$(info CUDA SDK version $(NVCC_VERSION) NOT SUPPORTED - (At least CUDA SDK 3.0 is required))

release_30: $(FOLDER_BUILD)/coreMyersGPU20.o $(FOLDER_BUILD)/myersGPU.o

release_50: $(FOLDER_BUILD)/coreMyersGPU35.o $(FOLDER_BUILD)/coreMyersGPU30.o \
            $(FOLDER_BUILD)/coreMyersGPU20.o $(FOLDER_BUILD)/myersGPU.o

release_60: $(FOLDER_BUILD)/coreMyersGPU50.o $(FOLDER_BUILD)/coreMyersGPU35.o \
            $(FOLDER_BUILD)/coreMyersGPU30.o $(FOLDER_BUILD)/coreMyersGPU20.o $(FOLDER_BUILD)/myersGPU.o
            
$(FOLDER_BUILD)/coreMyersGPU20.o: coreMyersGPU-CC20.cu
	nvcc $(NVCC_COMPILE_FLAGS) -gencode arch=compute_20,code=sm_20 -c coreMyersGPU-CC20.cu -o $(FOLDER_BUILD)/coreMyersGPU20.o

$(FOLDER_BUILD)/coreMyersGPU30.o: coreMyersGPU-CC30.cu
	nvcc $(NVCC_COMPILE_FLAGS) -gencode arch=compute_30,code=sm_30 -c coreMyersGPU-CC30.cu -o $(FOLDER_BUILD)/coreMyersGPU30.o

$(FOLDER_BUILD)/coreMyersGPU35.o: coreMyersGPU-CC35.cu
	nvcc $(NVCC_COMPILE_FLAGS) -gencode arch=compute_35,code=sm_35 -c coreMyersGPU-CC35.cu -o $(FOLDER_BUILD)/coreMyersGPU35.o

$(FOLDER_BUILD)/coreMyersGPU50.o: coreMyersGPU-CC50.cu
	nvcc $(NVCC_COMPILE_FLAGS) -gencode arch=compute_50,code=sm_50 -c coreMyersGPU-CC50.cu -o $(FOLDER_BUILD)/coreMyersGPU50.o
	
$(FOLDER_BUILD)/myersGPU.o: myersGPU.c
	gcc $(GCC_COMPILE_FLAGS) $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o $(FOLDER_BUILD)/myersGPU.o

##############################
#  MYERS TESTING VERSIONS	 #
##############################

profile: NVCC_COMPILE_FLAGS=-O3 -m64 --ptxas-options=-v -lineinfo
profile: GCC_COMPILE_FLAGS=-O3 -m64
profile: coreMyersGPU20 coreMyersGPU30 coreMyersGPU35 coreMyersGPU50
	gcc $(GCC_COMPILE_FLAGS) $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o

debug: NVCC_COMPILE_FLAGS=-G -O0 -m64 --ptxas-options=-v
debug: GCC_COMPILE_FLAGS=-g -O0 -m64
debug: myersGPU coreMyersGPU20 coreMyersGPU30 coreMyersGPU35 coreMyersGPU50

test: NVCC_COMPILE_FLAGS=-O3 -m64
test: GCC_COMPILE_FLAGS=-O3 -m64
test: myersGPU coreMyersGPU20 coreMyersGPU30 coreMyersGPU35 coreMyersGPU50
	gcc $(GCC_COMPILE_FLAGS) $(CUDA_LIBRARY_FLAGS) test.c $(OBJS) -o ../build/test -fopenmp

test-profile: NVCC_COMPILE_FLAGS=-O3 -m64 --ptxas-options=-v -lineinfo
test-profile: GCC_COMPILE_FLAGS=-O3 -m64
test-profile: myersGPU coreMyersGPU20 coreMyersGPU30 coreMyersGPU35 coreMyersGPU50
	gcc $(GCC_COMPILE_FLAGS) $(CUDA_LIBRARY_FLAGS) test.c $(OBJS) -o ../build/test-profile -fopenmp

test-debug: NVCC_COMPILE_FLAGS=-O0 -G -m64 --ptxas-options=-v -lineinfo
test-debug: GCC_COMPILE_FLAGS=-O0 -g -m64 --ptxas-options=-v -lineinfo
test-debug: myersGPU coreMyersGPU20 coreMyersGPU30 coreMyersGPU35 coreMyersGPU50
	gcc $(GCC_COMPILE_FLAGS) $(CUDA_LIBRARY_FLAGS) test.c $(OBJS) -o ../build/test-debug -fopenmp

coreMyersGPU20:
	nvcc $(NVCC_COMPILE_FLAGS) arch=compute_20,code=sm_20 -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	
coreMyersGPU30:
	nvcc $(NVCC_COMPILE_FLAGS) arch=compute_30,code=sm_30 -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	
coreMyersGPU35:
	nvcc $(NVCC_COMPILE_FLAGS) arch=compute_35,code=sm_35 -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	
coreMyersGPU50:
	nvcc $(NVCC_COMPILE_FLAGS) arch=compute_50,code=compute_50,sm_50 -c coreMyersGPU-CC50.cu -o ../build/coreMyersGPU50.o
	
myersGPU:
	gcc $(GCC_COMPILE_FLAGS) $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o

clean:
	rm -f ../build/*.o

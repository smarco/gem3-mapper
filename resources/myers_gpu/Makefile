SHELL := /bin/bash

NVCC_VERSION = $(shell nvcc --version | grep release | sed 's/.*release //' |  sed 's/,.*//' | sed 's/\.//')
CUDA_LIBRARY_FLAGS = -L/usr/local/cuda/lib64 -I/usr/local/cuda/include -lcuda -lcudart

NVCC_MINIMAL_VERSION = $(shell [[ ($(NVCC_VERSION) -ge "0")  && ($(NVCC_VERSION) -le  "29") ]] && echo "10"; \
				   		       [[ ($(NVCC_VERSION) -ge "30") && ($(NVCC_VERSION) -le  "49") ]] && echo "30"; \
			   				   [[ ($(NVCC_VERSION) -ge "50") && ($(NVCC_VERSION) -le  "59") ]] && echo "50"; \
			   				   [[ ($(NVCC_VERSION) -ge "60") && 						    ]] && echo "60";)

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

release_10:
	$(info CUDA SDK version $(NVCC_VERSION) NOT SUPPORTED - (At least CUDA SDK 3.0 is required))

release_30:
	nvcc -O3 -m64 -gencode=arch=compute_20,code=\’sm_20,compute_20\’ -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	gcc  -O3 -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o

release_50:
	nvcc -O3 -m64 -gencode=arch=compute_20,code=\sm_20,compute_20\ -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	nvcc -O3 -m64 -gencode=arch=compute_30,code=\sm_30,compute_30\ -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	nvcc -O3 -m64 -gencode=arch=compute_35,code=\sm_35,compute_35\ -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	gcc  -O3 -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o

release_60:
	nvcc -O3 -m64 -gencode arch=compute_20,code=\"sm_20,compute_20\" -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	nvcc -O3 -m64 -gencode arch=compute_30,code=\"sm_30,compute_30\" -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	nvcc -O3 -m64 -gencode arch=compute_35,code=\"sm_35,compute_35\" -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	nvcc -O3 -m64 -gencode arch=compute_50,code=\"sm_50,compute_50\" -c coreMyersGPU-CC50.cu -o ../build/coreMyersGPU50.o
	gcc  -O3 -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o


##############################
#  MYERS TESTING VERSIONS	 #
##############################

profile:
	nvcc -O3 -m64 arch=compute_20,code=sm_20 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	nvcc -O3 -m64 arch=compute_30,code=sm_30 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	nvcc -O3 -m64 arch=compute_35,code=sm_35 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	nvcc -O3 -m64 arch=compute_50,code=compute_50 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC50.cu -o ../build/coreMyersGPU50.o
	gcc  -O3 -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o

debug:
	nvcc -G -O0 -m64 arch=compute_20,code=sm_20 --ptxas-options=-v -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	nvcc -G -O0 -m64 arch=compute_30,code=sm_30 --ptxas-options=-v -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	nvcc -G -O0 -m64 arch=compute_35,code=sm_35 --ptxas-options=-v -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	nvcc -G -O0 -m64 arch=compute_50,code=compute_50,sm_50 --ptxas-options=-v -c coreMyersGPU-CC50.cu -o ../build/coreMyersGPU50.o
	gcc  -g -O0 -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o

test:
	nvcc -O3 -m64 -arch=compute_20 -code=sm_20 -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	nvcc -O3 -m64 -arch=compute_30 -code=sm_30 -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	nvcc -O3 -m64 -arch=compute_35 -code=sm_35 -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	nvcc -O3 -m64 -arch=compute_50 -code=compute_50,sm_50 -c coreMyersGPU-CC50.cu -o ../build/coreMyersGPU50.o
	gcc  -O3 -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o
	gcc  -O3 -m64 $(CUDA_LIBRARY_FLAGS) test.c ../build/myersGPU.o ../build/coreMyersGPU20.o ../build/coreMyersGPU30.o ../build/coreMyersGPU35.o ../build/coreMyersGPU50.o -o ../build/test -fopenmp

test-profile:
	nvcc -O3 -m64 -arch=compute_20 -code=sm_20 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	nvcc -O3 -m64 -arch=compute_30 -code=sm_30 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	nvcc -O3 -m64 -arch=compute_35 -code=sm_35 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	nvcc -O3 -m64 -arch=compute_50 -code=compute_50,sm_50 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC50.cu -o ../build/coreMyersGPU50.o
	gcc  -O3 -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o
	gcc  -O3 -m64 $(CUDA_LIBRARY_FLAGS) test.c ../build/myersGPU.o ../build/coreMyersGPU20.o ../build/coreMyersGPU30.o ../build/coreMyersGPU35.o ../build/coreMyersGPU50.o -o ../build/test-profile -fopenmp

test-debug:
	nvcc -O0 -G -m64 -arch=compute_20 -code=sm_20 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC20.cu -o ../build/coreMyersGPU20.o
	nvcc -O0 -G -m64 -arch=compute_30 -code=sm_30 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC30.cu -o ../build/coreMyersGPU30.o
	nvcc -O0 -G -m64 -arch=compute_35 -code=sm_35 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC35.cu -o ../build/coreMyersGPU35.o
	nvcc -O0 -G -m64 -arch=compute_50 -code=compute_50,sm_50 --ptxas-options=-v -lineinfo -c coreMyersGPU-CC50.cu -o ../build/coreMyersGPU50.o
	gcc  -O0 -g -m64 $(CUDA_LIBRARY_FLAGS) -c myersGPU.c -o ../build/myersGPU.o
	gcc  -O0 -g -m64 $(CUDA_LIBRARY_FLAGS) test.c ../build/myersGPU.o ../build/coreMyersGPU20.o ../build/coreMyersGPU30.o ../build/coreMyersGPU35.o ../build/coreMyersGPU50.o -o ../build/test-debug -fopenmp


clean:
	rm -f ../build/*.o

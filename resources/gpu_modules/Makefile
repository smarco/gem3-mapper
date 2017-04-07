#
#  GEM-Cutter "Highly optimized genomic resources for GPUs"
#  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
#
#  Licensed under GNU General Public License 3.0 or later.
#  Some rights reserved. See LICENSE, AUTHORS.
#  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
#
 
include ../../Makefile.mk

SHELL:=/bin/bash
CUDA_LIBRARY_FLAGS=$(CUDA_PATH_INCLUDE) $(CUDA_PATH_LIB) $(CUDA_LIB)
FOLDER_BUILD=build
FOLDER_SOURCE=src
FOLDER_BIN=bin
FOLDER_TOOLS=tools

NVCC_VERSION=$(shell $(NVCC) --version | grep release | sed 's/.*release //' |  sed 's/,.*//' | sed 's/\.//')
NVCC_REQUIRED_VERSION = $(shell [[ ($(NVCC_VERSION) -lt "50") ]] && echo "false";)
NVCC_MINIMAL_VERSION  = $(shell [[ ($(NVCC_VERSION) -ge "0")  && ($(NVCC_VERSION) -lt  "30") ]] && echo "10"; \
                                [[ ($(NVCC_VERSION) -ge "30") && ($(NVCC_VERSION) -lt  "50") ]] && echo "20"; \
                                [[ ($(NVCC_VERSION) -ge "50") && ($(NVCC_VERSION) -lt  "65") ]] && echo "30"; \
                                [[ ($(NVCC_VERSION) -ge "65") && ($(NVCC_VERSION) -lt  "70") ]] && echo "50"; \
                                [[ ($(NVCC_VERSION) -ge "70") && ($(NVCC_VERSION) -lt  "80") ]] && echo "52"; \
                                [[ ($(NVCC_VERSION) -ge "80") ]] && echo "60";)
ifeq ($(NVCC_REQUIRED_VERSION), false)
$(error CUDA SDK version $(NVCC_VERSION) NOT SUPPORTED - (At least CUDA SDK 5.0 is required))
endif

$(shell test -d $(FOLDER_BUILD) || mkdir $(FOLDER_BUILD))
$(shell test -d $(FOLDER_BIN) || mkdir $(FOLDER_BIN))

## Force to test JIT compiler:
## export CUDA_FORCE_PTX_JIT=1
CUDA_SASS_FLAG_20=-gencode arch=compute_20,code=\"sm_20,compute_20\" -gencode arch=compute_20,code=sm_21
CUDA_SASS_FLAG_30=$(CUDA_SASS_FLAG_20) -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\"
CUDA_SASS_FLAG_37=$(CUDA_SASS_FLAG_30) -gencode arch=compute_37,code=sm_37
CUDA_SASS_FLAG_50=$(CUDA_SASS_FLAG_37) -gencode arch=compute_50,code=\"sm_50,compute_50\"
CUDA_SASS_FLAG_52=$(CUDA_SASS_FLAG_50) -gencode arch=compute_52,code=\"sm_52,compute_52\"
CUDA_SASS_FLAG_60=$(CUDA_SASS_FLAG_52) -gencode arch=compute_60,code=\"sm_60,compute_60\" -gencode arch=compute_61,code=\"sm_61,compute_61\" -gencode arch=compute_62,code=\"sm_62,compute_62\"
CUDA_SASS_FLAGS=$(CUDA_SASS_FLAG_$(NVCC_MINIMAL_VERSION))

CUDA_MODULES=gpu_fmi_decode gpu_fmi_ssearch gpu_fmi_asearch gpu_bpm_align gpu_bpm_filter gpu_kmer_filter gpu_sa_decode
CUDA_SRCS=$(addprefix $(FOLDER_SOURCE)/, $(addsuffix .cu, $(CUDA_MODULES)))
CUDA_OBJS=$(addprefix $(FOLDER_BUILD)/, $(addsuffix .o, $(CUDA_MODULES)))

BASICS=gpu_commons gpu_buffer gpu_errors gpu_io gpu_sample gpu_module gpu_devices gpu_index gpu_reference
FMI_MODULES=gpu_fmi_index gpu_fmi_table gpu_fmi_primitives gpu_fmi_primitives_decode gpu_fmi_primitives_ssearch gpu_fmi_primitives_asearch
SA_MODULES=gpu_sa_index gpu_sa_primitives
BPM_MODULES=gpu_bpm_primitives_filter gpu_bpm_primitives_align
KMER_MODULES=gpu_kmer_primitives_filter
MODULES= $(FMI_MODULES) $(SA_MODULES) $(BPM_MODULES) $(KMER_MODULES) $(BASICS)
SRCS=$(addprefix $(FOLDER_SOURCE)/, $(addsuffix .c, $(MODULES)))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(addsuffix .o, $(MODULES)))

TOOLS=gpu_benchmark_bpm gpu_benchmark_decode gpu_benchmark_search
TOOLS_SRC=$(addprefix $(FOLDER_TOOLS)/, $(addsuffix .c, $(TOOLS)))
TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(TOOLS))

BUILDERS=gpu_build_index gpu_build_reference gpu_build_regions
BUILDERS_SRC=$(addprefix $(FOLDER_TOOLS)/, $(addsuffix .c, $(BUILDERS)))
BUILDERS_BIN=$(addprefix $(FOLDER_BIN)/, $(BUILDERS))

all: release

release: NVCC_COMPILE_FLAGS=-O3 -m64 -Xptxas="-dlcm=ca"
release: GCC_COMPILE_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL)
release: $(OBJS) $(CUDA_OBJS) 

devel: NVCC_COMPILE_FLAGS=-g -O3 -m64 -lineinfo -Xptxas="-dlcm=ca"
devel: GCC_COMPILE_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_DEVEL)
devel: $(OBJS) $(CUDA_OBJS)

profile: NVCC_COMPILE_FLAGS=-O3 -m64 --ptxas-options=-v -lineinfo -Xptxas="-dlcm=ca"
profile: GCC_COMPILE_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_DEVEL) $(FLAGS_PROFILE)
profile: $(OBJS) $(CUDA_OBJS) 
	
debug: NVCC_COMPILE_FLAGS=-g -G -O0 -m64 -lineinfo -Xptxas="-dlcm=ca"
debug: GCC_COMPILE_FLAGS=$(FLAGS_DEBUG) $(FLAGS_GENERAL) $(FLAGS_DEVEL) $(FLAGS_PROFILE)
debug: $(OBJS) $(CUDA_OBJS) 
	
tools: release $(TOOLS_BIN) $(BUILDERS_BIN)
tools_profile: profile $(TOOLS_BIN) $(BUILDERS_BIN)
tools_devel: devel $(TOOLS_BIN) $(BUILDERS_BIN)
tools_debug: debug $(TOOLS_BIN) $(BUILDERS_BIN)


$(FOLDER_BUILD)/%.o: $(FOLDER_SOURCE)/%.cu
	$(NVCC) $(NVCC_COMPILE_FLAGS) $(CUDA_SASS_FLAGS) -c $< -o $@
	
$(FOLDER_BUILD)/%.o: $(FOLDER_SOURCE)/%.c
	$(CC) $(GCC_COMPILE_FLAGS) -c $< -o $@ $(CUDA_LIBRARY_FLAGS) -lrt

link: $(OBJS) $(CUDA_OBJS)
	ld -r $(OBJS) $(CUDA_OBJS) -o $(FOLDER_BUILD)/gem_gpu.o

$(FOLDER_BIN)/gpu_benchmark_%: $(FOLDER_TOOLS)/gpu_benchmark_%.c
	$(CC) $(GCC_COMPILE_FLAGS) $(FOLDER_BUILD)/*.o $< -o $@ $(CUDA_LIBRARY_FLAGS) -fopenmp -lrt

$(FOLDER_BIN)/gpu_build_%: $(FOLDER_TOOLS)/gpu_build_%.c
	$(CC) $(GCC_COMPILE_FLAGS) $< -o $@ -lrt

clean:
	rm -f $(FOLDER_BIN)/* $(FOLDER_BUILD)/*.o

#============================================================================
# PROJECT: GEM library
# FILE: Makefile
# DATE: 04/09/2014
# DESCRIPTION: 
#============================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

all: release

release: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=myers_gpu
endif
debug: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=myers_gpu all
endif
setup: 
	@mkdir -p build

clean:
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=myers_gpu clean
endif
	@rm -rf build

#============================================================================
# PROJECT: GEM library
# FILE: Makefile
# DATE: 04/09/2014
# DESCRIPTION: 
#============================================================================

# Definitions
ROOT_PATH=..
include $(ROOT_PATH)/Makefile.mk

all: devel

release: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gpu_modules release
endif

devel: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gpu_modules devel
endif

profile: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gpu_modules profile
endif

debug: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gpu_modules debug
endif

setup: 
	@mkdir -p build

clean:
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gpu_modules clean
endif
	@rm -rf build

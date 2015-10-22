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
	$(MAKE) --directory=myers_gpu release
endif

devel: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=myers_gpu release
endif

profile: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=myers_gpu profile
endif

debug: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=myers_gpu debug
endif

setup: 
	@mkdir -p build

clean:
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=myers_gpu clean
endif
	@rm -rf build

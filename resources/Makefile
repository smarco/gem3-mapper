#============================================================================
# PROJECT: GEM library
# FILE: Makefile
# DATE: 04/09/2014
# DESCRIPTION: 
#============================================================================

# Definitions
ROOT_PATH=..
GEM_CUTTER=$(if $(shell ls -A gem-cutter),1,0)
include $(ROOT_PATH)/Makefile.mk

all: devel

release: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gem-cutter release
	cp -u gem-cutter/build/*.o ../build
endif

devel: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gem-cutter devel
	cp -u gem-cutter/build/*.o ../build
endif

profile: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gem-cutter profile
	cp -u gem-cutter/build/*.o ../build
endif

debug: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gem-cutter debug
	cp -u gem-cutter/build/*.o ../build
endif

clean: setup
ifeq ($(and $(HAVE_CUDA),$(GEM_CUTTER)),)
	$(MAKE) --directory=gem-cutter clean
endif
	@rm -rf build

setup: 
	@mkdir -p build
ifeq ($(GEM_CUTTER),0)
	@git submodule update --init --recursive
	GEM_CUTTER=1
endif


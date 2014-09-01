#============================================================================
# PROJECT: GEM library
# FILE: Makefile
# DATE: 31/10/2012
# DESCRIPTION: 
#============================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

all: release

release: setup
	$(MAKE) --directory=myers_gpu

debug: setup
  $(MAKE) --directory=myers_gpu all
  
setup: 
	@mkdir -p build

clean:
	$(MAKE) --directory=myers_gpu clean
	@rm -rf build

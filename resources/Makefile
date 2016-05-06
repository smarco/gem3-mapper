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
	cp -u gpu_modules/build/*.o ../build
endif

devel: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gpu_modules devel
	cp -u gpu_modules/build/*.o ../build
endif

profile: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gpu_modules profile
	cp -u gpu_modules/build/*.o ../build
endif

debug: setup
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gpu_modules debug
	cp -u gpu_modules/build/*.o ../build
endif

setup: 
	@mkdir -p build

clean:
ifeq ($(HAVE_CUDA),1)
	$(MAKE) --directory=gpu_modules clean
endif
	@rm -rf build

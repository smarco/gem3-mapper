#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 31/10/2012
# DESCRIPTION: Top-level makefile
#============================================================================

# Definitions
ROOT_PATH=$(CURDIR)
include Makefile.mk

all: release

release: setup
	$(MAKE) --directory=resources
	$(MAKE) --directory=src
	$(MAKE) --directory=tools
	
debug: setup
	$(MAKE) --directory=resources
	$(MAKE) --directory=src debug
	$(MAKE) --directory=tools debug
	
static:	setup
  $(MAKE) --directory=resources
	$(MAKE) --directory=src
	$(MAKE) --directory=tools static
	
profile: setup
  $(MAKE) --directory=resources
	$(MAKE) --directory=src profile
	$(MAKE) --directory=tools profile

check: setup debug
	$(MAKE) --directory=test check
	
setup: 
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)

clean:
	$(MAKE) --directory=resources clean
	$(MAKE) --directory=test clean
	@rm -rf $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)

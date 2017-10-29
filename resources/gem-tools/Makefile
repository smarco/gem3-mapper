#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Top level makefile
#============================================================================

# Definitions
ROOT_PATH=$(CURDIR)
include Makefile.mk

all: release

static:	setup
	$(MAKE) --directory=src
	$(MAKE) --directory=resources/src
	$(MAKE) --directory=tools static

release: setup	
	$(MAKE) --directory=src
	$(MAKE) --directory=resources/src
	$(MAKE) --directory=tools

debug: setup
	$(MAKE) --directory=src debug
	$(MAKE) --directory=resources/src debug
	$(MAKE) --directory=tools debug

check: setup debug
	$(MAKE) --directory=test check
	
setup: 
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)

clean:
	$(MAKE) --directory=test clean
	@rm -rf $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)
	

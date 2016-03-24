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

# Optimized 
release: setup
	$(MAKE) --directory=resources release
	$(MAKE) --directory=src release
	
# Optimized + Static // TODO
static:	setup
	$(MAKE) --directory=resources release
	$(MAKE) --directory=src release
	
# Optimized + DebugSymbols
devel: setup
	$(MAKE) --directory=resources devel
	$(MAKE) --directory=src devel
	
# Optimized + DebugSymbols + GEMProfile
profile: setup
	$(MAKE) --directory=resources profile
	$(MAKE) --directory=src profile
	
# DebugSymbols + GEMProfile + GEMDebug
debug: setup
	$(MAKE) --directory=resources debug
	$(MAKE) --directory=src debug
	
complexity:
	pmccabe -vt include/*/*.h src/*/*.c
	
setup: 
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)
	@ln -s $(FOLDER_RESOURCES) $(FOLDER_INCLUDE)/resources 2> /dev/null ||:

clean:
	$(MAKE) --directory=resources clean
	@rm -rf $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)
	@rm -rf $(FOLDER_INCLUDE)/resources

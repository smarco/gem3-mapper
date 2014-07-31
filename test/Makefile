#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# DESCRIPTION: Builds and runs all the tests
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

FOLDER_TEST_BUILD=./build
FOLDER_TEST_SCRIPT=./script

UTESTS=utest_commons
ITESTS=

UTESTS_FLAGS=$(ARCH_FLAGS) $(DEBUG_FLAGS)
COVERAGE_FLAGS=-g -Wall -fprofile-arcs -ftest-coverage $(TESTS_FLAGS)

LIBS=-lgemcore_c -lpthread -lcheck -lm $(LIBS_ZLIB) $(LIBS_BZLIB)

all: check coverage
check: clean setup $(UTESTS) $(ITESTS) end_banner
coverage: clean setup $(COVERAGE) end_banner

$(UTESTS):
	$(CC) $(UTESTS_FLAGS) $(INCLUDE_FLAGS) $(LIB_PATH_FLAGS) -o $(FOLDER_TEST_BUILD)/$@ $@.c $(LIBS)
	@echo "=======================================================================>>"
	@echo "==>> Testing " $@
	@echo "=======================================================================>>"
	-$(FOLDER_TEST_BUILD)/$@
	
$(ITESTS): 
	/bin/bash $(FOLDER_TEST_SCRIPT)/$@.sh
	
$(COVERAGE): clean setup
	$(CC) $(COVERAGE_FLAGS) $(INCLUDE_FLAGS) $(LIB_PATH_FLAGS) -o $(FOLDER_TEST_BUILD)/$@ $@.c $(LIBS)
	$(FOLDER_TEST_BUILD)/$@
	
setup:
	@mkdir -p $(FOLDER_TEST_BUILD) $(FOLDER_TEST_SCRIPT)
	
end_banner: 	
	@echo "=======================================================================>>"
	
clean:
	@rm -rf $(FOLDER_TEST_BUILD)

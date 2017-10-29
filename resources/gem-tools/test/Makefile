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
FOLDER_TEST_REPORTS=./reports

GT_UTESTS=gt_utest_commons gt_utest_core_structures gt_utest_parsers
GT_ITESTS=gt_itest_map_parser

GT_UTESTS_FLAGS=$(ARCH_FLAGS) $(DEBUG_FLAGS)
GT_COVERAGE_FLAGS=-g -Wall -fprofile-arcs -ftest-coverage $(GT_TESTS_FLAGS)

LIBS=-lpthread -lgemtools -lcheck -lz -lbz2 -fopenmp

all: check coverage

check: clean setup $(GT_UTESTS) $(GT_ITESTS) end_banner

coverage: clean setup $(GT_COVERAGE) end_banner

$(GT_UTESTS):
	$(CC) $(GT_UTESTS_FLAGS) $(INCLUDE_FLAGS) $(LIB_PATH_FLAGS) -o $(FOLDER_TEST_BUILD)/$@ $@.c $(LIBS)
	@echo "=======================================================================>>"
	@echo "==>> Testing " $@
	@echo "=======================================================================>>"
	-$(FOLDER_TEST_BUILD)/$@
	
$(GT_ITESTS): 
	/bin/bash $(FOLDER_TEST_SCRIPT)/$@.sh
	
$(GT_COVERAGE): clean setup
	$(CC) $(GT_COVERAGE_FLAGS) $(INCLUDE_FLAGS) $(LIB_PATH_FLAGS) -o $(FOLDER_TEST_BUILD)/$@ $@.c $(LIBS)
	$(FOLDER_TEST_BUILD)/$@
	
setup:
	@mkdir -p $(FOLDER_TEST_BUILD) $(FOLDER_TEST_SCRIPT) $(FOLDER_TEST_REPORTS)
	
end_banner: 	
	@echo "=======================================================================>>"
	
clean:
	@rm -rf $(FOLDER_TEST_BUILD) $(FOLDER_TEST_REPORTS)

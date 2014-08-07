#============================================================================
# PROJECT: GEMMapper
# FILE: Makefile
# DATE: 31/10/2012
# DESCRIPTION:
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

GEM_TOOLS=gem-mapper # gem-constructor gem-indexer 
GEM_TOOLS_SRC=$(addsuffix .c, $(GEM_TOOLS))
GEM_TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(GEM_TOOLS))
LIBS=-lgemcore_c -lpthread -lm $(LIBS_ZLIB) $(LIBS_BZLIB)

all: release

release: GEM_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(SUPPRESS_CHECKS) 
release: $(GEM_TOOLS_BIN)

static: GEM_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(SUPPRESS_CHECKS) -static
static: $(GEM_TOOLS_BIN)

profile: GEM_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(SUPPRESS_CHECKS) -g -ggdb3 
profile: LIBS+=/opt/intel/vtune_amplifier_xe_2013/lib64/libittnotify.a -ldl
profile: $(GEM_TOOLS_BIN)

debug: GEM_FLAGS=-O0 $(GENERAL_FLAGS) $(DEBUG_FLAGS)
debug: $(GEM_TOOLS_BIN)

$(FOLDER_BIN)/gem-mapper:
	$(CC) $(GEM_FLAGS) -o $@ $(notdir $@).c $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS) $(CUDA_LIBRARY_FLAGS) -fopenmp

# $(GEM_TOOLS_BIN): $(FOLDER_LIB)/libgemcore_c.a $(GEM_TOOLS_SRC)
#	$(CC) $(GEM_FLAGS) -o $@ $(notdir $@).c $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS)
	
	
	
	

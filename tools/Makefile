#============================================================================
# PROJECT: GEMMapper
# FILE: Makefile
# DATE: 31/10/2012
# DESCRIPTION:
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

GEM_TOOLS=gem-mapper gem-indexer gem-constructor
GEM_TOOLS_SRC=$(addsuffix .c, $(GEM_TOOLS))
GEM_TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(GEM_TOOLS))
LIBS=-lgemcore_c -lpthread -lm $(LIBS_ZLIB) $(LIBS_BZLIB)
ifeq ($(OS),Linux)
LIBS+=-lrt
endif

all: devel

release: GEM_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(SUPPRESS_CHECKS) -DGEM_NOPROFILE
release: $(GEM_TOOLS_BIN)

static: GEM_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(SUPPRESS_CHECKS) -static
static: $(GEM_TOOLS_BIN)

devel: GEM_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS)
devel: $(GEM_TOOLS_BIN)

profile: GEM_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS)
profile: LIBS+=/opt/intel/vtune_amplifier_xe_2013/lib64/libittnotify.a -ldl
profile: $(GEM_TOOLS_BIN)

debug: GEM_FLAGS=-O0 $(GENERAL_FLAGS) $(DEBUG_FLAGS) $(GEM_DEBUG_FLAGS)
debug: $(GEM_TOOLS_BIN)

$(FOLDER_BIN)/gem-mapper: $(GEM_CORE_LIB) gem-mapper.c
	$(CC) $(GEM_FLAGS) $(notdir $@).c -o $@ $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS) $(DEF_CUDA) $(LIBS_CUDA) 

$(FOLDER_BIN)/gem-indexer: $(GEM_CORE_LIB) gem-indexer.c
	$(CC) $(GEM_FLAGS) $(notdir $@).c -o $@ $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS)

$(FOLDER_BIN)/gem-constructor: $(GEM_CORE_LIB) gem-constructor.c
	$(CC) $(GEM_FLAGS) $(notdir $@).c -o $@ $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS)
	
$(FOLDER_BIN)/gem-filtering: $(GEM_CORE_LIB) gem-filtering.c
	$(CC) $(GEM_FLAGS) $(notdir $@).c -o $@ $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS)
	

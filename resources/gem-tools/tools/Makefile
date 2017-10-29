#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds the GT-tools
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

GEM_TOOLS=gt.stats gt.filter gt.mapset gt.map2sam gt.checkCandidates gt.constructor gt.pe2se
GEM_TOOLS_SRC=$(addsuffix .c, $(GEM_TOOLS))
GEM_TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(GEM_TOOLS))
LIBS=-lgemtools -lpthread -lm $(LIBS_OPENMP) $(LIBS_ZLIB) $(LIBS_BZLIB)

GT_FLAGS=$(GENERAL_FLAGS) $(DEF_BZLIB) $(DEF_ZLIB) $(DEF_OPENMP) $(ARCH_FLAGS) 
GT_OPT_FLAGS=$(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)

all: GEM_TOOLS_FLAGS=-O4 $(GT_FLAGS) $(GT_OPT_FLAGS) 
all: $(GEM_TOOLS_BIN)

static: GEM_TOOLS_FLAGS=-O4 $(GT_FLAGS) $(GT_OPT_FLAGS) -static
static: $(GEM_TOOLS_BIN)

debug: GEM_TOOLS_FLAGS=-O0 $(GT_FLAGS) $(DEBUG_FLAGS)
debug: $(GEM_TOOLS_BIN)

$(FOLDER_BIN)/gt.stats: LIBS:=$(LIBS) -ljson
$(GEM_TOOLS_BIN): $(FOLDER_LIB)/libgemtools.a $(GEM_TOOLS_SRC)
	$(CC) $(GEM_TOOLS_FLAGS) -o $@ $(notdir $@).c $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS)

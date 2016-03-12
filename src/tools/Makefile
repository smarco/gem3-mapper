#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 31/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds all the sources
#==================================================================================================

###############################################################################
# Definitions
###############################################################################
ROOT_PATH=../..
include $(ROOT_PATH)/Makefile.mk

###############################################################################
# Tools
###############################################################################
GEM_TOOLS=gem-mapper gem-indexer gem-constructor gem-retriever
GEM_TOOLS_SRC=$(addsuffix .c, $(GEM_TOOLS))
GEM_TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(GEM_TOOLS))

###############################################################################
# Compilation Modes
###############################################################################
all: devel

release: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_LINK) $(FLAGS_SUPPRESS_CHECKS)
release: $(GEM_TOOLS_BIN)

static: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_LINK) $(FLAGS_SUPPRESS_CHECKS) -static
static: $(GEM_TOOLS_BIN)

devel: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_LINK) $(FLAGS_DEBUG)
devel: LIBS+=$(VTUNE_PROFILE_LIB)
devel: $(GEM_TOOLS_BIN)

profile: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_LINK) $(FLAGS_DEBUG) $(FLAGS_GEM_DEBUG) $(FLAGS_PROFILE)
profile: LIBS+=$(CUDA_PROFILE_LIB)
profile: LIBS+=$(VTUNE_PROFILE_LIB)
profile: PATH_INCLUDE+=$(CUDA_PATH_INCLUDE)
profile: $(GEM_TOOLS_BIN)

debug: GEM_FLAGS=-O0 $(FLAGS_GENERAL) $(FLAGS_DEBUG) $(FLAGS_GEM_DEBUG) $(FLAGS_PROFILE)
debug: LIBS+=$(CUDA_PROFILE_LIB)
debug: PATH_INCLUDE+=$(CUDA_PATH_INCLUDE)
debug: $(GEM_TOOLS_BIN)

###############################################################################
# Rules
###############################################################################

# GEM-library
$(LIB_GEM_CORE): $(FOLDER_BUILD)/*.o
ifeq ($(HAVE_CUDA),1)
	$(AR) -rcs $(LIB_GEM_CORE) $(FOLDER_BUILD)/*.o $(FOLDER_GEMGPU_BUILD)/*.o
else
	$(AR) -rcs $(LIB_GEM_CORE) $(FOLDER_BUILD)/*.o
endif
	
# GEM Tools	
$(FOLDER_BIN)/gem-mapper: $(LIB_GEM_CORE) gem-mapper.c 
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB) $(CUDA_LIB) $(LIBS) 

$(FOLDER_BIN)/gem-indexer: $(LIB_GEM_CORE) gem-indexer.c
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB) $(CUDA_LIB) $(LIBS) 
	
$(FOLDER_BIN)/gem-retriever: $(LIB_GEM_CORE) gem-retriever.c
	$(CC) $(GEM_FLAGS) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB) $(CUDA_LIB) $(LIBS) 

$(FOLDER_BIN)/gem-constructor: $(LIB_GEM_CORE) gem-constructor.c
	$(CC) $(GEM_FLAGS) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB) $(CUDA_LIB) $(LIBS) 
	


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
# Version
###############################################################################
VERSION=$(shell if [ -d $(ROOT_PATH)/.git ]; then git describe --abbrev=4 --dirty --always; else cat $(ROOT_PATH)/VERSION; fi)

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

release: GEM_VERSION=-DGEM_VERSION=$(VERSION)-release 
release: GEM_PATHS=$(PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB)
release: GEM_LIBS=$(CUDA_LIB) $(LIBS) $(VTUNE_PROFILE_LIB)
release: $(GEM_TOOLS_BIN)

devel: GEM_VERSION=-DGEM_VERSION=$(VERSION)-devel
devel: GEM_PATHS=$(PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB)
devel: GEM_LIBS=$(CUDA_LIB) $(LIBS) $(VTUNE_PROFILE_LIB)
devel: $(GEM_TOOLS_BIN)

profile: GEM_VERSION=-DGEM_VERSION=$(VERSION)-profile
profile: GEM_PATHS=$(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB)
profile: GEM_LIBS=$(CUDA_LIB) $(LIBS) $(CUDA_PROFILE_LIB) $(VTUNE_PROFILE_LIB)
profile: $(GEM_TOOLS_BIN)

debug: GEM_VERSION=-DGEM_VERSION=$(VERSION)-debug 
debug: GEM_PATHS=$(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB) 
debug: GEM_LIBS=$(CUDA_LIB) $(LIBS) $(CUDA_PROFILE_LIB) $(VTUNE_PROFILE_LIB)
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
	$(CC) $(GEM_FLAGS) $(GEM_VERSION) $(DEF_CUDA) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(FLAGS_LINK) $(GEM_PATHS) $(GEM_LIBS)  

$(FOLDER_BIN)/gem-indexer: $(LIB_GEM_CORE) gem-indexer.c
	$(CC) $(GEM_FLAGS) $(GEM_VERSION) $(DEF_CUDA) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(FLAGS_LINK) $(GEM_PATHS) $(GEM_LIBS)  
	
$(FOLDER_BIN)/gem-retriever: $(LIB_GEM_CORE) gem-retriever.c
	$(CC) $(GEM_FLAGS) $(GEM_VERSION) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(FLAGS_LINK) $(GEM_PATHS) $(GEM_LIBS) 

$(FOLDER_BIN)/gem-constructor: $(LIB_GEM_CORE) gem-constructor.c
	$(CC) $(GEM_FLAGS) $(GEM_VERSION) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(FLAGS_LINK) $(GEM_PATHS) $(GEM_LIBS) 
	


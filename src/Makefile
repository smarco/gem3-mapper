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
ROOT_PATH=..
include $(ROOT_PATH)/Makefile.mk



###############################################################################
# Compile all Modules
###############################################################################
SUBDIRS=align \
        approximate_search \
        archive \
        data_structures \
        filtering \
        fm_index \
        gpu \
        io \
        mapper \
        matches \
        neighborhood_search \
        search_pipeline \
        stats \
        system \
        utils \
        tools

release: MODE=release
release: $(SUBDIRS)

devel: MODE=devel
devel: $(SUBDIRS)

profile: MODE=profile
profile: $(SUBDIRS)

debug: MODE=debug
debug: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) --directory=$@ $(MODE)

###############################################################################
# GEM-library
###############################################################################

$(LIB_GEM_CORE): $(OBJS)
ifeq ($(HAVE_CUDA),1)
	$(AR) -rcs $(LIB_GEM_CORE) $(FOLDER_BUILD)/*.o $(FOLDER_GEMGPU_BUILD)/*.o
else
	$(AR) -rcs $(LIB_GEM_CORE) $(FOLDER_BUILD)/*.o
endif


.PHONY: $(SUBDIRS)

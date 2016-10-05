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
        filtering \
        fm_index \
        gpu \
        io \
        mapper \
        matches \
        neighborhood_search \
        profiler \
        search_pipeline \
        stats \
        system \
        text \
        utils \
        tools

release: GEM_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_SUPPRESS_CHECKS)
release: MODE=release
release: $(SUBDIRS)

devel: GEM_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_DEVEL) 
devel: MODE=devel
devel: $(SUBDIRS)

profile: GEM_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_DEVEL) $(FLAGS_PROFILE) $(FLAGS_GEM_DEBUG) $(FLAGS_GEM_PROFILE) 
profile: MODE=profile
profile: $(SUBDIRS)

debug: GEM_FLAGS=$(FLAGS_DEBUG) $(FLAGS_GENERAL) $(FLAGS_DEVEL) $(FLAGS_PROFILE) $(FLAGS_GEM_DEBUG) $(FLAGS_GEM_PROFILE)
debug: MODE=debug
debug: $(SUBDIRS)

# Subdir rule
export GEM_FLAGS
export MODE
$(SUBDIRS):
	$(MAKE) --directory=$@ $(MODE)

.PHONY: $(SUBDIRS)

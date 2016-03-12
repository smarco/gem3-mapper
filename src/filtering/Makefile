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
# Modules
###############################################################################
MODULES=filtering_candidates_align \
        filtering_candidates_align_local \
        filtering_candidates \
        filtering_candidates_extend \
        filtering_candidates_process_buffered \
        filtering_candidates_process \
        filtering_candidates_verify_buffered \
        filtering_candidates_verify \
        filtering_region_align \
        filtering_region_align_configure \
        filtering_region \
        filtering_region_cache \
        filtering_region_verify \
        kmer_counting \
        region_profile_adaptive \
        region_profile \
        region_profile_fixed \
        region_profile_schedule \
        region_set

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

###############################################################################
# Compilation Modes
###############################################################################
all: devel

release: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_SUPPRESS_CHECKS)
release: $(OBJS)

devel: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_DEBUG)
devel: $(OBJS)

profile: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_DEBUG) $(FLAGS_GEM_DEBUG) $(FLAGS_PROFILE)
profile: $(OBJS)

debug: GEM_FLAGS=-O0 $(FLAGS_GENERAL) $(FLAGS_DEBUG) $(FLAGS_GEM_DEBUG) $(FLAGS_PROFILE)
debug: $(OBJS)

###############################################################################
# Rules
###############################################################################
	
# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_FLAGS) $(PATH_INCLUDE) -c $< -o $@

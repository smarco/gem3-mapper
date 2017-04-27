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
ROOT_PATH=../../..
include $(ROOT_PATH)/Makefile.mk

###############################################################################
# Modules
###############################################################################
MODULES=filtering_candidates \
        filtering_candidates_accessors \
        filtering_candidates_align \
        filtering_candidates_align_local \
        filtering_candidates_buffered \
        filtering_candidates_buffered_bpm_align \
        filtering_candidates_buffered_bpm_distance \
        filtering_candidates_buffered_kmer_filter \
        filtering_candidates_buffered_process \
        filtering_candidates_classify \
        filtering_candidates_extend \
        filtering_candidates_process \
        filtering_candidates_verify

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

###############################################################################
# Compilation Modes
###############################################################################
all: devel

release: $(OBJS)

devel: $(OBJS)

profile: $(OBJS)

debug: $(OBJS)

###############################################################################
# Rules
###############################################################################
	
# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_FLAGS) $(PATH_INCLUDE) -c $< -o $@

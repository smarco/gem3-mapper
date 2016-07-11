#==============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 31/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds all the sources
#==============================================================================

###############################################################################
# Definitions
###############################################################################
ROOT_PATH=../..
include $(ROOT_PATH)/Makefile.mk

###############################################################################
# Modules
###############################################################################
MODULES=match_align \
        match_alignment \
        match_align_normalize \
        match_align_rl \
        match_align_swg \
        match_align_swg_local \
        matches \
        matches_cigar \
        matches_classify \
        matches_classify_logit \
        matches_classify_logit_models \
        matches_counters \
        matches_metrics \
        matches_predictors \
        match_scaffold \
        match_scaffold_compose \
        match_scaffold_levenshtein \
        match_scaffold_region_chain \
        paired_map \
        paired_matches

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

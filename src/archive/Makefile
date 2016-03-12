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
MODULES=archive_builder \
        archive_builder_index \
        archive_builder_text \
        archive_builder_text_parser \
        archive \
        archive_check \
        archive_score_pe \
        archive_score_se \
        archive_search \
        archive_search_cache \
        archive_search_paired_parameters \
        archive_search_parameters \
        archive_search_pe \
        archive_search_pe_stepwise \
        archive_search_se \
        archive_search_se_stepwise \
        archive_select \
        archive_select_parameters \
        archive_text \
        archive_text_rl \
        locator_builder \
        locator \
        sampled_rl

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

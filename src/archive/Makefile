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
        archive_search_pe_stages \
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

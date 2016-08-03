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
MODULES=approximate_search_control \
        approximate_search_filtering_adaptive \
        approximate_search_filtering_complete \
        approximate_search_generate_candidates \
        approximate_search_hybrid \
        approximate_search_neighborhood \
        approximate_search_region_profile \
        approximate_search_stages \
        approximate_search_stepwise \
        approximate_search_verify_candidates \
        approximate_search \
        mm_search

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

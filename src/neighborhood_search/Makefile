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
MODULES=dp_matrix \
        nsearch_hamming \
        nsearch_levenshtein \
        nsearch_levenshtein_scheduled \
        nsearch_levenshtein_state \
        nsearch_partition \
        nsearch_operation \
        nsearch_schedule

SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

###############################################################################
# Compilation Modes
###############################################################################
all: devel

release: GEM_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_SUPPRESS_CHECKS)
release: $(OBJS)

devel: GEM_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_DEVEL)
devel: $(OBJS)

profile: GEM_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_DEVEL) $(FLAGS_PROFILE) $(FLAGS_GEM_DEBUG) $(FLAGS_GEM_PROFILE)
profile: $(OBJS)

debug: GEM_FLAGS=$(FLAGS_DEBUG) $(FLAGS_GENERAL) $(FLAGS_DEVEL) $(FLAGS_PROFILE) $(FLAGS_GEM_DEBUG) $(FLAGS_GEM_PROFILE)
debug: $(OBJS)

###############################################################################
# Rules
###############################################################################
	
# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_FLAGS) $(PATH_INCLUDE) -c $< -o $@

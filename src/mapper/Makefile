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
MODULES=mapper_bisulfite \
        mapper \
        mapper_cuda \
        mapper_cuda_pe \
        mapper_cuda_se \
        mapper_profile \
        mapper_profile_cuda \
        mapper_profile_misc \
        mapper_profile_search \
        mapper_stats

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

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
MODULES=profiler \
        profiler_cuda \
        profiler_gem \
        profiler_timer \
        profiler_vtune

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

# Building Exceptions
$(FOLDER_BUILD)/profiler_cuda.o : profiler_cuda.c
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(DEF_CUDA_PROFILE) $(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) -c $< -o $@
$(FOLDER_BUILD)/profiler_vtune.o : profiler_vtune.c
	$(CC) $(GEM_FLAGS) $(PATH_INCLUDE) $(FLAGS_GEM_VTUNE) $(VTUNE_PROFILE_INCLUDE) -c $< -o $@

# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_FLAGS) $(PATH_INCLUDE) -c $< -o $@

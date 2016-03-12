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
MODULES=commons \
        errors \
        fm \
        gruntime \
        gthread \
        mm \
        mm_pool \
        mm_slab \
        mm_stack \
        profiler \
        profiler_cuda \
        profiler_gem \
        profiler_timer \
        report

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
profile: DEF_CUDA_PROFILE=-DGEM_PROFILE_CUDA
profile: $(OBJS)

debug: GEM_FLAGS=-O0 $(FLAGS_GENERAL) $(FLAGS_DEBUG) $(FLAGS_GEM_DEBUG) $(FLAGS_PROFILE)
debug: $(OBJS)

###############################################################################
# Rules
###############################################################################

# Building Exceptions
$(FOLDER_BUILD)/fm.o : fm.c
	$(CC) $(GEM_FLAGS) $(DEF_BZLIB) $(DEF_ZLIB) $(PATH_INCLUDE) -c $< -o $@
$(FOLDER_BUILD)/profiler_cuda.o : profiler_cuda.c
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(DEF_CUDA_PROFILE) $(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) -c $< -o $@

# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_FLAGS) $(PATH_INCLUDE) -c $< -o $@

#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds all the sources of the GT-library, and the library itself (.a)
#==================================================================================================

# Definitions
ROOT_PATH=../..
include ../../Makefile.mk

MODULES=json
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))
RESOURCE_LIBS=$(FOLDER_LIB)/libjson.a

all: GEM_TOOLS_FLAGS=$(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(RESOURCE_LIBS)

debug: GEM_TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(RESOURCE_LIBS)

$(RESOURCE_LIBS): $(OBJS)
	$(AR) -rcs $(RESOURCE_LIBS) $(OBJS)

ifeq ($(HAVE_OPENMP),1)
OPENMP_FLAGS:= -fopenmp
else
OPENMP_FLAGS:=
endif

$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_TOOLS_FLAGS) $(INCLUDE_FLAGS) -c $< -o $@

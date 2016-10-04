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
SUBDIRS=builder \
        search \
        score

MODULES=archive \
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
all: release

release: $(SUBDIRS)
release: $(OBJS)

devel: $(SUBDIRS)
devel: $(OBJS)

profile: $(SUBDIRS)
profile: $(OBJS)

debug: $(SUBDIRS)
debug: $(OBJS)

###############################################################################
# Rules
###############################################################################

# Building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_FLAGS) $(PATH_INCLUDE) -c $< -o $@

# Subdir rule
export GEM_FLAGS
export MODE
$(SUBDIRS):
	$(MAKE) --directory=$@ $(MODE)
	
.PHONY: $(SUBDIRS)

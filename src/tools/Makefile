###############################################################################
#  GEM-Mapper v3 (GEM3)
#  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
#
#  This file is part of GEM-Mapper v3 (GEM3).
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# PROJECT: GEM-Mapper v3 (GEM3)
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
###############################################################################

###############################################################################
# Definitions
###############################################################################
ROOT_PATH=../..
include $(ROOT_PATH)/Makefile.mk

###############################################################################
# Version
###############################################################################
VERSION=$(shell if [ -d $(ROOT_PATH)/.git ]; then git describe --abbrev=4 --dirty --always; else cat $(ROOT_PATH)/VERSION; fi)

###############################################################################
# Tools
###############################################################################
GEM_TOOLS=gem-mapper gem-indexer gem-retriever
GEM_TOOLS_SRC=$(addsuffix .c, $(GEM_TOOLS))
GEM_TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(GEM_TOOLS))

###############################################################################
# Modules
###############################################################################
SUBDIRS=interface

###############################################################################
# Compilation Modes
###############################################################################
all: devel

release: $(SUBDIRS)
release: GEM_VERSION=-DGEM_VERSION=$(VERSION)-release 
release: GEM_PATHS=$(PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB)
release: GEM_LIBS=$(CUDA_LIB) $(LIBS) $(VTUNE_PROFILE_LIB)
release: $(GEM_TOOLS_BIN)

devel: $(SUBDIRS)
devel: GEM_VERSION=-DGEM_VERSION=$(VERSION)-devel
devel: GEM_PATHS=$(PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB)
devel: GEM_LIBS=$(CUDA_LIB) $(LIBS) $(VTUNE_PROFILE_LIB)
devel: $(GEM_TOOLS_BIN)

profile: $(SUBDIRS)
profile: GEM_VERSION=-DGEM_VERSION=$(VERSION)-profile
profile: GEM_PATHS=$(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB)
profile: GEM_LIBS=$(CUDA_LIB) $(LIBS) $(CUDA_PROFILE_LIB) $(VTUNE_PROFILE_LIB)
profile: $(GEM_TOOLS_BIN)

debug: $(SUBDIRS)
debug: GEM_VERSION=-DGEM_VERSION=$(VERSION)-debug 
debug: GEM_PATHS=$(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) $(PATH_LIB) $(CUDA_PATH_LIB) 
debug: GEM_LIBS=$(CUDA_LIB) $(LIBS) $(CUDA_PROFILE_LIB) $(VTUNE_PROFILE_LIB)
debug: $(GEM_TOOLS_BIN)

###############################################################################
# Rules
###############################################################################

# GEM-library
$(LIB_GEM_CORE): $(FOLDER_BUILD)/*.o
	$(AR) -rcs $(LIB_GEM_CORE) $(FOLDER_BUILD)/*.o
	
# GEM Tools	
$(FOLDER_BIN)/gem-mapper: $(FOLDER_BUILD)/*.o gem-mapper.c 
	$(CC) $(GEM_FLAGS) $(GEM_VERSION) $(DEF_CUDA) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(FLAGS_LINK) $(GEM_PATHS) $(GEM_LIBS)  

$(FOLDER_BIN)/gem-indexer: $(FOLDER_BUILD)/*.o gem-indexer.c
	$(CC) $(GEM_FLAGS) $(GEM_VERSION) $(DEF_CUDA) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(FLAGS_LINK) $(GEM_PATHS) $(GEM_LIBS)
	
$(FOLDER_BIN)/gem-mappability: $(FOLDER_BUILD)/*.o gem-mappability.c 
	$(CC) $(GEM_FLAGS) $(GEM_VERSION) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(FLAGS_LINK) $(GEM_PATHS) $(GEM_LIBS)  
	
$(FOLDER_BIN)/gem-retriever: $(FOLDER_BUILD)/*.o gem-retriever.c
	$(CC) $(GEM_FLAGS) $(GEM_VERSION) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(FLAGS_LINK) $(GEM_PATHS) $(GEM_LIBS) 

$(FOLDER_BIN)/gem-constructor: $(FOLDER_BUILD)/*.o gem-constructor.c
	$(CC) $(GEM_FLAGS) $(GEM_VERSION) $(notdir $@).c $(FOLDER_BUILD)/*.o -o $@ $(FLAGS_LINK) $(GEM_PATHS) $(GEM_LIBS) 
	
# Subdir rule
export GEM_FLAGS
export MODE
$(SUBDIRS):
	$(MAKE) --directory=$@ $(MODE)
	
.PHONY: $(SUBDIRS)

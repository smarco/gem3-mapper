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
ROOT_PATH=..
include $(ROOT_PATH)/Makefile.mk

###############################################################################
# Compile all Modules
###############################################################################
SUBDIRS=align \
        approximate_search \
        archive \
        filtering \
        fm_index \
        gpu \
        io \
        mapper \
        matches \
        neighborhood_search \
        profiler \
        search_pipeline \
        stats \
        system \
        text \
        utils \
        tools

release: GEM_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_SUPPRESS_CHECKS)
release: MODE=release
release: $(SUBDIRS)

devel: GEM_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_DEVEL) 
devel: MODE=devel
devel: $(SUBDIRS)

profile: GEM_FLAGS=$(FLAGS_OPT) $(FLAGS_GENERAL) $(FLAGS_DEVEL) $(FLAGS_PROFILE) $(FLAGS_GEM_DEBUG) $(FLAGS_GEM_PROFILE) 
profile: MODE=profile
profile: $(SUBDIRS)

debug: GEM_FLAGS=$(FLAGS_DEBUG) $(FLAGS_GENERAL) $(FLAGS_DEVEL) $(FLAGS_PROFILE) $(FLAGS_GEM_DEBUG) $(FLAGS_GEM_PROFILE)
debug: MODE=debug
debug: $(SUBDIRS)

# Subdir rule
export GEM_FLAGS
export MODE
$(SUBDIRS): 
	$(MAKE) --directory=$@ $(MODE)

.PHONY: $(SUBDIRS)

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

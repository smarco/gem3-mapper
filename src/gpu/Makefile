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
MODULES=gpu_buffer_bpm_align \
        gpu_buffer_bpm_distance \
        gpu_buffer_bpm_pattern \
        gpu_buffer_collection \
        gpu_buffer_fmi_ssearch \
        gpu_buffer_fmi_asearch \
        gpu_buffer_fmi_decode \
        gpu_buffer_kmer_filter \
        gpu_config \
        gpu_structures

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
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) -c $< -o $@

#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 31/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds all the sources of the GEMTools and the library (.a)
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

# Basic Modules
MODULES =commons gthread profiler report errors gthread \
         fm mm mm_slab mm_stack mm_pool \
         vector ihash shash priority_queue string_buffer segmented_vector \
         profiler stats_vector stats_matrix
# I/O
MODULES+=input_file input_buffer buffered_input_file input_parser \
         input_fasta_parser input_multifasta_parser input_graph_parser \
         output_file output_buffer buffered_output_file \
         output_map output_sam
# Advance Data Structures
MODULES+=sequence dna_text cdna_text cdna_bitwise_text \
         sparse_bitmap sparse_array_locator packed_integer_array \
         text_collection graph_text graph_text_builder \
# GEM Structures
MODULES+=archive_builder archive_text_builder archive_hypertext_builder archive_index_builder sa_builder \
         bwt rank_mtable sampled_sa sampled_rl fm_index interval_set matches locator \
         archive archive_search archive_search_group mm_search \
         archive_select archive_select_parameters archive_text_retrieve 
# GEM Search Modules
MODULES+=approximate_search_parameters approximate_search neighborhood_search matches_align \
         filtering_candidates region_profile bpm_align bpm_align_gpu quality_model
# GEM Top-level Modules
MODULES+=mapper mapper_cuda gem_core
# Misc
MODULES+=mapper_profile options_menu 
        
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))
GEM_CORE_LIB=$(FOLDER_LIB)/libgemcore_c.a

all: release

release: COMPILE_C_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(SUPPRESS_CHECKS)
release: $(GEM_CORE_LIB)

profile: COMPILE_C_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(SUPPRESS_CHECKS) -g -ggdb3 -I/opt/intel/vtune_amplifier_xe_2013/include/
profile: $(GEM_CORE_LIB)

debug: COMPILE_C_FLAGS=-O0 $(GENERAL_FLAGS) $(DEBUG_FLAGS)
debug: $(GEM_CORE_LIB)

$(GEM_CORE_LIB): $(OBJS)
ifeq ($(HAVE_CUDA),1)
	$(AR) -rcs $(GEM_CORE_LIB) $(FOLDER_BUILD)/*.o $(FOLDER_RESOURCES_BUILD)/*.o
else
	$(AR) -rcs $(GEM_CORE_LIB) $(FOLDER_BUILD)/*.o
endif

# Building Exceptions (Dependencies)
$(FOLDER_BUILD)/fm.o : fm.c
	$(CC) $(COMPILE_C_FLAGS) $(DEF_BZLIB) $(DEF_ZLIB) $(INCLUDE_FLAGS) -c $< -o $@
$(FOLDER_BUILD)/bpm_align_gpu.o : bpm_align_gpu.c
	$(CC) $(COMPILE_C_FLAGS) $(DEF_CUDA) $(INCLUDE_FLAGS) -c $< -o $@
	
# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(COMPILE_C_FLAGS) $(INCLUDE_FLAGS) -c $< -o $@

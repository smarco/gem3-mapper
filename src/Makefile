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
MODULES=commons gthread report errors
# Basic Data Structures
MODULES+=vector ihash shash priority_queue string_buffer segmented_vector
# Stats Data Structures
MODULES+=profiler stats_vector stats_matrix profiler_cuda
# Memory/Files Modules
MODULES+=fm mm mm_slab mm_stack mm_pool
# Input Data Structures
MODULES+=input_file input_buffer buffered_input_file
# Output Data Structures
MODULES+=output_file output_buffer buffered_output_file
# BM Data Structures
MODULES+=sparse_bitmap sparse_array_locator packed_integer_array
# Sequence Data Structures
MODULES+=sequence dna_text cdna_text cdna_bitwise_text text_collection
# BWT Structures
MODULES+=bwt_commons bwt_sbasic bwt_basic
# FM-Index Structures
MODULES+=fm_index fm_index_search fm_index_query \
         locator locator_builder \
         rank_mtable rank_mtable_builder  \
         sampled_sa sampled_rl sa_builder \
         interval_set region_set
# GEM Matches
MODULES+=matches matches_counters matches_metrics matches_classify 
# GEM paired-Matches
MODULES+=paired_matches paired_map paired_matches_classify
# GEM Pattern Alignment
MODULES+=pattern bpm_align swg_align bpm_align_gpu match_scaffold match_align
# GEM Region Profile
MODULES+=region_profile region_profile_adaptive \
         region_profile_boost region_profile_schedule \
         region_profile_fixed
# GEM Filtering
MODULES+=filtering_candidates filtering_region kmer_counting 
# GEM Search Support Modules
MODULES+=search_parameters paired_search_parameters select_parameters quality_model mm_search
# GEM Approximate Search (Filtering)
MODULES+=approximate_search approximate_search_filtering_base \
         approximate_search_filtering_stages approximate_search_filtering_control \
         approximate_search_filtering_adaptive approximate_search_neighborhood \
         approximate_search_filtering_complete
# GEM Approximate Search (Neighborhood Search)
MODULES+=nsearch nsearch_hamming nsearch_levenshtein nsearch_partition nsearch_schedule
# GEM Archive
MODULES+=archive archive_builder archive_builder_index \
         archive_builder_text archive_builder_text_parser \
         archive_text archive_search archive_search_group \
         archive_search_se archive_search_pe \
         archive_select archive_score
# Input Parsers
MODULES+=input_parser input_fasta_parser input_multifasta_parser
# Output Printers
MODULES+=output_sam output_map
# GEM Support Modules
MODULES+=options_menu
# Report Stats Module
MODULES+=report_stats
# GEM Main Modules
MODULES+=mapper mapper_cuda mapper_bisulfite mapper_stats mapper_profile gem_core
        
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

all: devel

release: COMPILE_C_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(SUPPRESS_CHECKS)
release: $(GEM_CORE_LIB)

devel: COMPILE_C_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS)
devel: $(GEM_CORE_LIB)

profile: COMPILE_C_FLAGS=$(GENERAL_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS) -I/opt/intel/vtune_amplifier_xe_2013/include/
profile: $(GEM_CORE_LIB)

debug: COMPILE_C_FLAGS=-O0 $(GENERAL_FLAGS) $(DEBUG_FLAGS) $(GEM_DEBUG_FLAGS)
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

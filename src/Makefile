#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 31/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds all the sources of the GEMTools and the library (.a)
#==================================================================================================

###############################################################################
# Definitions
###############################################################################
ROOT_PATH=..
include $(ROOT_PATH)/Makefile.mk

###############################################################################
# Modules
###############################################################################
# Basic Modules
MODULES=commons gthread report errors
# Basic Data Structures
MODULES+=vector ihash shash priority_queue string_buffer segmented_vector dp_matrix
# Stats Data Structures
MODULES+=stats_vector stats_matrix stats_counter
# Profile Data Structures
MODULES+=profiler profiler_timer profiler_gem profiler_cuda
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
MODULES+=matches matches_cigar \
         matches_counters matches_metrics \
         matches_classify 
# GEM Paired-Matches
MODULES+=paired_matches paired_map paired_matches_classify
# GEM Pattern Alignment
MODULES+=pattern align \
         align_bpm align_bpm_pattern \
         align_bpm_distance align_ond \
         align_swg align_swg_score \
         align_swg_simd match_scaffold
# GEM Match Alignment
MODULES+=match_alignment match_align \
         match_align_swg match_align_local
# GEM Region Profile
MODULES+=region_profile region_profile_adaptive \
         region_profile_boost region_profile_schedule \
         region_profile_fixed
# GEM Filtering
MODULES+=filtering_candidates filtering_candidates_verify \
         filtering_candidates_align filtering_candidates_extend \
         filtering_candidates_process filtering_candidates_process_buffered \
         filtering_region_verify filtering_candidates_verify_buffered \
         filtering_region filtering_region_align \
         filtering_region_cache kmer_counting 
# GEM Approximate Search (Filtering)
MODULES+=approximate_search approximate_search_filtering_adaptive \
         approximate_search_filtering_complete approximate_search_neighborhood \
         approximate_search_filtering_stages approximate_search_control \
         approximate_search_region_profile approximate_search_generate_candidates \
         approximate_search_verify_candidates approximate_search_stepwise
# GEM Approximate Search (Neighborhood Search)
MODULES+=nsearch nsearch_partition \
         nsearch_schedule nsearch_hamming \
         nsearch_levenshtein nsearch_levenshtein_state
# GEM Archive Builder
MODULES+=archive archive_builder archive_builder_index \
         archive_builder_text archive_builder_text_parser
# GEM Archive Search Support
MODULES+=archive_text archive_search_cache \
         archive_select archive_select_parameters \
         archive_score archive_check \
         mm_search quality_model \
# GEM Archive Search
MODULES+=archive_search \
         archive_search_se archive_search_parameters \
         archive_search_pe archive_search_paired_parameters \
         archive_search_se_stepwise archive_search_pe_stepwise
# GEM GPU-Buffers
MODULES+=gpu_buffer_collection gpu_buffer_align_bpm \
         gpu_buffer_fmi_bsearch gpu_buffer_fmi_decode
# GEM Search-Pipeline & Stages
MODULES+=search_pipeline \
         search_stage_region_profile search_stage_region_profile_buffer \
         search_stage_decode_candidates search_stage_decode_candidates_buffer \
         search_stage_verify_candidates search_stage_verify_candidates_buffer
# Input Parsers
MODULES+=input_parser input_fasta_parser input_multifasta_parser
# Output Printers
MODULES+=output_sam output_map
# GEM Support Modules
MODULES+=options_menu
# Stats Modules
MODULES+=mapper_stats mapper_profile \
         mapper_profile_cuda mapper_profile_search \
         mapper_profile_misc report_stats
# GEM Main Modules
MODULES+=mapper mapper_cuda mapper_cuda_se mapper_cuda_pe mapper_bisulfite gem_core
        
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

###############################################################################
# Compilation Modes
###############################################################################
all: devel

release: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_SUPPRESS_CHECKS)
release: $(LIB_GEM_CORE)

devel: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_DEBUG)
devel: $(LIB_GEM_CORE)

profile: GEM_FLAGS=$(FLAGS_GENERAL) $(FLAGS_OPT) $(FLAGS_DEBUG) $(FLAGS_PROFILE)
profile: DEF_CUDA_PROFILE=-DGEM_PROFILE_CUDA
profile: $(LIB_GEM_CORE)

debug: GEM_FLAGS=-O0 $(FLAGS_GENERAL) $(FLAGS_DEBUG) $(FLAGS_PROFILE)
debug: $(LIB_GEM_CORE)

###############################################################################
# Rules
###############################################################################
# GEM-library
$(LIB_GEM_CORE): $(OBJS)
ifeq ($(HAVE_CUDA),1)
	$(AR) -rcs $(LIB_GEM_CORE) $(FOLDER_BUILD)/*.o $(FOLDER_GEMGPU_BUILD)/*.o
else
	$(AR) -rcs $(LIB_GEM_CORE) $(FOLDER_BUILD)/*.o
endif

# Building Exceptions (Dependencies)
$(FOLDER_BUILD)/fm.o : fm.c
	$(CC) $(GEM_FLAGS) $(DEF_BZLIB) $(DEF_ZLIB) $(PATH_INCLUDE) -c $< -o $@
$(FOLDER_BUILD)/gpu_buffer_collection.o : gpu_buffer_collection.c
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) -c $< -o $@
$(FOLDER_BUILD)/gpu_buffer_fmi_bsearch.o : gpu_buffer_fmi_bsearch.c
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) -c $< -o $@
$(FOLDER_BUILD)/gpu_buffer_fmi_decode.o : gpu_buffer_fmi_decode.c
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) -c $< -o $@
$(FOLDER_BUILD)/gpu_buffer_align_bpm.o : gpu_buffer_align_bpm.c
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) -c $< -o $@
$(FOLDER_BUILD)/profiler_cuda.o : profiler_cuda.c
	$(CC) $(GEM_FLAGS) $(DEF_CUDA) $(DEF_CUDA_PROFILE) $(PATH_INCLUDE) $(CUDA_PATH_INCLUDE) -c $< -o $@
	
# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_FLAGS) $(PATH_INCLUDE) -c $< -o $@

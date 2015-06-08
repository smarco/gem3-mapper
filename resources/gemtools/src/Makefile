#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds all the sources of the GT-library, and the library itself (.a)
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

MODULES=gem_tools \
        gt_commons gt_error gt_mm gt_fm gt_profiler \
        gt_ihash gt_shash gt_vector gt_string gt_json \
        gt_attributes gt_dna_string gt_dna_read gt_compact_dna_string \
        gt_template gt_alignment gt_map gt_misms \
        gt_template_utils gt_alignment_utils gt_counters_utils \
        gt_map_metrics gt_map_score gt_map_utils \
        gt_map_align gt_map_align_bpm gt_map_align_bpm_simd \
        gt_sequence_archive gt_segmented_sequence gt_gemIdx_loader \
        gt_input_file gt_buffered_input_file \
        gt_input_parser gt_input_map_parser gt_input_fasta_parser gt_input_generic_parser \
        gt_input_map_utils \
        gt_input_sam_parser gt_sam_attributes \
        gt_buffered_output_file gt_output_file gt_generic_printer gt_output_buffer \
        gt_output_printer gt_output_map gt_output_fasta gt_output_sam gt_output_generic_printer \
        gt_stats gt_options_menu
        
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))
GT_LIB=$(FOLDER_LIB)/libgemtools.a

all: GEM_TOOLS_FLAGS=$(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(GT_LIB)

debug: GEM_TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(GT_LIB)

$(GT_LIB): $(OBJS)
	$(AR) -rcs $(GT_LIB) $(FOLDER_BUILD)/*.o

$(FOLDER_BUILD)/gt_fm.o : gt_fm.c
	$(CC) $(GEM_TOOLS_FLAGS) $(DEF_BZLIB) $(DEF_ZLIB) $(INCLUDE_FLAGS) -c $< -o $@
$(FOLDER_BUILD)/gt_options_menu.o : gt_options_menu.c
	$(CC) $(GEM_TOOLS_FLAGS) $(DEF_OPENMP) $(INCLUDE_FLAGS) -c $< -o $@
$(FOLDER_BUILD)/gt_map_align_bpm_simd.o : gt_map_align_bpm_simd.c
	$(CC) $(GEM_TOOLS_FLAGS) -msse4.2 $(DEF_OPENMP) $(INCLUDE_FLAGS) -c $< -o $@
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_TOOLS_FLAGS) $(INCLUDE_FLAGS) -c $< -o $@

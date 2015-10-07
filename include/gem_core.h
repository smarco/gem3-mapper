/*
 * PROJECT: GEMMapper
 * FILE: gem_core.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef GEM_CORE_H_
#define GEM_CORE_H_

/*
 * Version
 */
#define GEM_CORE_VERSION 3.1
#define GEM_CORE_VERSION_STR "3.1"

// GEM essentials
#include "essentials.h"
#include "segmented_vector.h"

// GEM Structures
#include "sparse_array_locator.h"
#include "packed_integer_array.h"
#include "rank_mtable.h"

#include "cdna_text.h"
#include "cdna_bitwise_text.h"
#include "dna_text.h"

// GEM Index
#include "sa_builder.h"
#include "bwt.h"
#include "fm_index.h"
#include "locator.h"
#include "archive.h"
#include "archive_builder.h"
#include "archive_builder_text.h"
#include "archive_builder_text_parser.h"
#include "archive_builder_index.h"
#include "select_parameters.h"
#include "archive_select.h"

// I/O
#include "input_file.h"
#include "input_parser.h"
#include "input_fasta_parser.h"
#include "output_map.h"
#include "output_sam.h"

// Profile
#include "profiler.h"
#include "mapper_profile.h"

// Options Menu (Adaptors + Helpers)
#include "options_menu.h"

// Mapper
#include "mapper.h"
#include "mapper_cuda.h"

// Report Stats
#include "report_stats.h"

/*
 * GEM Runtime
 */
void gem_runtime_init(const uint64_t num_threads,const uint64_t max_memory,char* const tmp_folder);
void gem_runtime_destroy();

#endif /* GEM_CORE_H_ */

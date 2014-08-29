/*
 * PROJECT: GEMMapper
 * FILE: gem_core.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef GEM_CORE_H_
#define GEM_CORE_H_

// GEM essentials
#include "essentials.h"
#include "segmented_vector.h"

// GEM Structures
#include "sparse_array_locator.h"
#include "packed_integer_array.h"
#include "rank_mtable.h"

#include "dna_string.h"
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

// I/O
#include "input_file.h"
#include "input_parser.h"
#include "input_fasta_parser.h"
#include "output_map.h"
#include "output_sam.h"

// Stats
#include "stats_vector.h"

// Options Menu (Adaptors + Helpers)
#include "options_menu.h"

// Mapper
#include "mapper.h"
#include "mapper_cuda.h"

/*
 * GEM Runtime
 */
GEM_INLINE void gem_runtime_init(const uint64_t max_memory,char* const tmp_folder,report_function_t report_function);

GEM_INLINE void gem_runtime_delete();

#endif /* GEM_CORE_H_ */

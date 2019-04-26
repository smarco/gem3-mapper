/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Mapper module encapsulates and provides accessors to all
 *   the parameters used by the mapper
 */

#include "mapper/mapper_parameters.h"

/*
 * Mapper Parameters
 */
void mapper_parameters_set_defaults_io(mapper_parameters_io_t* const io) {
  const uint64_t num_processors = system_get_num_processors();
  /* Input */
  io->index_file_name=NULL;
  io->separated_input_files=false;
  io->input_file_name=NULL;
  io->input_file_name_end1=NULL;
  io->input_file_name_end2=NULL;
  io->input_compression=FM_REGULAR_FILE;
  /* I/O */
  io->input_block_size = BUFFER_SIZE_32M;
  io->input_num_blocks = num_processors;
  io->input_buffer_size = BUFFER_SIZE_4M;
  io->output_file_name=NULL;
  io->output_compression=FM_REGULAR_FILE;
  /* I/O Attributes */
  io->fastq_strictly_normalized = false;
  /* Output */
  io->output_format = SAM;
  output_sam_parameters_set_defaults(&io->sam_parameters);
  output_map_parameters_set_defaults(&io->map_parameters);
  io->output_buffer_size = BUFFER_SIZE_4M;
  io->output_num_buffers = 10*num_processors; // Lazy allocation
  io->report_file_name = NULL;
  io->mapper_ticker_step = 100000;
}
void mapper_parameters_set_defaults_system(mapper_parameters_system_t* const system) {
  /* System */
  const uint64_t num_processors = system_get_num_processors();
  system->num_threads=num_processors;
  system->max_memory=0;
  system->tmp_folder=NULL;
}
void mapper_parameters_set_defaults_cuda(mapper_parameters_cuda_t* const cuda) {
  /* CUDA settings */
  const uint64_t num_processors = system_get_num_processors();
  /* CUDA */
  cuda->gpu_enabled=false;
  cuda->gpu_devices=UINT64_MAX;
  /* I/O */
  cuda->input_block_size = BUFFER_SIZE_32M;
  cuda->input_buffer_size = BUFFER_SIZE_4M;
  cuda->output_buffer_size = BUFFER_SIZE_4M;
  cuda->output_num_buffers = 10*num_processors; // Lazy allocation
  /* BPM Buffers */
  cuda->gpu_buffer_size = CONVERT_B_TO_MB(BUFFER_SIZE_1M);
  cuda->num_fmi_bsearch_buffers = 2;
  cuda->num_fmi_decode_buffers = 3;
  cuda->num_kmer_filter_buffers = 3;
  cuda->num_bpm_distance_buffers = 3;
  cuda->num_bpm_align_buffers = 3;
  /* Stages Configuration */
  cuda->cpu_emulation=false;
}
//void mapper_parameters_set_defaults_hints(mapper_parameters_hints_t* const hints) {
//  /* Hints */
//}
void mapper_parameters_set_defaults_misc(mapper_parameters_misc_t* const misc) {
  /* QC */
  misc->quality_control = false;
  misc->profile = false;
  misc->profile_reduce_type = reduce_sample;
  /* Verbose */
  misc->verbose_user=true;
  misc->verbose_dev=false;
}
void mapper_parameters_set_defaults(mapper_parameters_t* const mapper_parameters) {
  /* CMD line */
  mapper_parameters->argc = 0;
  mapper_parameters->argv = NULL;
  /* GEM Structures */
  mapper_parameters->archive = NULL;
  mapper_parameters->input_file = NULL;
  mapper_parameters->input_file_end1 = NULL;
  mapper_parameters->input_file_end2 = NULL;
  MUTEX_INIT(mapper_parameters->input_file_mutex);
  mapper_parameters->output_stream = NULL;
  mapper_parameters->output_file = NULL;
  /* Mapper Type */
  mapper_parameters->mapper_type = mapper_se;
  /* Stats Report */
  mapper_parameters->global_mapping_stats = NULL;
  /* I/O Parameters */
  mapper_parameters_set_defaults_io(&mapper_parameters->io);
  /* Search Parameters (single-end/paired-end) */
  search_parameters_init(&mapper_parameters->search_parameters);
  /* System */
  mapper_parameters_set_defaults_system(&mapper_parameters->system);
  /* CUDA settings */
  mapper_parameters_set_defaults_cuda(&mapper_parameters->cuda);
  /* Miscellaneous */
  mapper_parameters_set_defaults_misc(&mapper_parameters->misc);
  /* Parsing Arguments Flags */
  mapper_parameters->min_reported_strata_set = false;
  mapper_parameters->max_reported_matches_set = false;
  mapper_parameters->bs_suffix1 = NULL;
  mapper_parameters->bs_suffix2 = NULL;
}
/*
 * Mapper parameters display
 */
void mapper_parameters_print(FILE* const stream,mapper_parameters_t* const parameters) {
  tab_fprintf(stream,"[GEM]>Mapper.parameters\n");
  /* CMD line */
  uint64_t i;
  tab_fprintf(stream,"  => Application %s\n",parameters->argv[0]);
  tab_fprintf(stream,"  => Arguments   ");
  for (i=1;i<parameters->argc;++i) {
    fprintf(stream,"%s ",parameters->argv[i]);
  }
  fprintf(stream,"\n");
}

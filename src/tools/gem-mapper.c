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
 */

#define REPORT_STATS

#include "utils/essentials.h"
#include "stats/report_stats.h"
#include "mapper/mapper.h"
#include "mapper/mapper_parameters.h"
#include "mapper/mapper_io.h"
#include "mapper/mapper_cuda.h"
#include "mapper/mapper_profile.h"
#include "mapper/mapper_profile_cuda.h"
#include "tools/interface/mapper_arguments.h"

/*
 * Version
 */
#define GEM_VERSION_STRING(version) QUOTE(version)
char* const gem_version = GEM_VERSION_STRING(GEM_VERSION);

/*
 * GEM-mapper I/O related functions
 */
input_file_sliced_t* gem_mapper_open_input_file(
    char* const input_file_name,const fm_type input_compression,
    const uint64_t input_block_size,const uint64_t input_num_blocks,
	const bool verbose_user) {
  // Open input file
  if(input_file_name != NULL) {
	  size_t l = strlen(input_file_name);
	  if(l > 1) {
		  if(input_file_name[l - 1] == '|') {
			  char *tmp = mm_malloc(l);
			  memcpy(tmp, input_file_name, l - 1);
			  tmp[l - 1] = 0;
			  gem_cond_log(verbose_user,"[Opening pipe from '%s']",tmp);
			  input_file_sliced_t *fp = input_file_sliced_popen(tmp,input_num_blocks,input_block_size);
			  mm_free(tmp);
			  return fp;
		  } else {
			  gem_cond_log(verbose_user,"[Opening input file '%s']",input_file_name);
			  return input_file_sliced_open(input_file_name,input_num_blocks,input_block_size);
		  }
	  }
  }
  gem_cond_log(verbose_user,"[Reading input file from stdin]");
	switch (input_compression) {
	 case FM_GZIPPED_FILE:
		return input_gzip_stream_sliced_open(stdin,input_num_blocks,input_block_size);
	 case FM_BZIPPED_FILE:
		return input_bzip_stream_sliced_open(stdin,input_num_blocks,input_block_size);
	 default:
		return input_stream_sliced_open(stdin,input_num_blocks,input_block_size);
  }

}
void gem_mapper_open_input(mapper_parameters_t* const parameters) {
  if (parameters->io.separated_input_files) {
    parameters->input_file_end1 = gem_mapper_open_input_file(
        parameters->io.input_file_name_end1,parameters->io.input_compression,
        parameters->io.input_block_size,parameters->io.input_num_blocks,
        parameters->misc.verbose_user);
    parameters->input_file_end2 = gem_mapper_open_input_file(
        parameters->io.input_file_name_end2,parameters->io.input_compression,
        parameters->io.input_block_size,parameters->io.input_num_blocks,
        parameters->misc.verbose_user);
  } else {
    parameters->input_file = gem_mapper_open_input_file(
        parameters->io.input_file_name,parameters->io.input_compression,
        parameters->io.input_block_size,parameters->io.input_num_blocks,
        parameters->misc.verbose_user);
  }
}
void gem_mapper_close_input(mapper_parameters_t* const parameters) {
  if (parameters->io.separated_input_files) {
    input_file_sliced_close(parameters->input_file_end1);
    input_file_sliced_close(parameters->input_file_end2);
  } else {
    input_file_sliced_close(parameters->input_file);
  }
}
void gem_mapper_open_output(mapper_parameters_t* const parameters) {
  // Open output stream
  if (parameters->io.output_file_name==NULL) {
    gem_cond_log(parameters->misc.verbose_user,"[Outputting to stdout]");
    parameters->output_stream = stdout;
  } else {
    gem_cond_log(parameters->misc.verbose_user,"[Outputting to '%s']",parameters->io.output_file_name);
    parameters->output_stream = fopen(parameters->io.output_file_name,"w");
    mapper_cond_error_msg(parameters->output_stream==NULL,
        "Couldn't open output file '%s'",parameters->io.output_file_name);
  }
  // Open output file
  const mapper_parameters_cuda_t* const cuda = &parameters->cuda;
  const uint64_t max_output_buffers = (cuda->gpu_enabled) ? cuda->output_num_buffers : parameters->io.output_num_buffers;
  const uint64_t output_buffer_size = (cuda->gpu_enabled) ? cuda->output_buffer_size : parameters->io.output_buffer_size;
  switch (parameters->io.output_compression) {
    case FM_GZIPPED_FILE:
      parameters->output_file = output_gzip_stream_new(parameters->output_stream,max_output_buffers,output_buffer_size);
      break;
    case FM_BZIPPED_FILE:
      parameters->output_file = output_bzip_stream_new(parameters->output_stream,max_output_buffers,output_buffer_size);
      break;
    default:
      parameters->output_file = output_stream_new(parameters->output_stream,max_output_buffers,output_buffer_size);
      break;
  }
}
void gem_mapper_close_output(mapper_parameters_t* const parameters) {
  output_file_close(parameters->output_file);
}
void gem_mapper_print_profile(mapper_parameters_t* const parameters) {
  // Reduce Stats
  switch (parameters->misc.profile_reduce_type) {
    case reduce_sum: PROF_REDUCE_SUM(); break;
    case reduce_max: PROF_REDUCE_MAX(); break;
    case reduce_min: PROF_REDUCE_MIN(); break;
    case reduce_mean: PROF_REDUCE_MEAN(); break;
    case reduce_sample: PROF_REDUCE_SAMPLE(); break;
    default: GEM_INVALID_CASE(); break;
  }
  // System
  system_print_info(gem_log_get_stream());
  // Parameters
  // mapper_parameters_print(gem_log_get_stream(),parameters);
  // Mapper
  if (!parameters->cuda.gpu_enabled) {
    // CPU Mapper
    switch (parameters->mapper_type) {
      case mapper_se:
        mapper_profile_print_mapper_se(gem_log_get_stream(),
            parameters->io.output_format==MAP,parameters->system.num_threads);
        break;
      case mapper_pe:
        mapper_profile_print_mapper_pe(gem_log_get_stream(),
            parameters->io.output_format==MAP,parameters->system.num_threads);
        break;
      default:
        break;
    }
  } else {
    // CUDA Mapper
    switch (parameters->mapper_type) {
      case mapper_se:
        mapper_profile_print_mapper_se_cuda(gem_log_get_stream(),
            parameters->io.output_format==MAP,parameters->system.num_threads);
        break;
      case mapper_pe:
        mapper_profile_print_mapper_pe_cuda(gem_log_get_stream(),
            parameters->io.output_format==MAP,parameters->system.num_threads);
        break;
      default:
        break;
    }
  }
}
/*
 * Main
 */
int main(int argc,char** argv) {
  // Parsing command-line options
  mapper_parameters_t parameters;
  mapper_parameters_set_defaults(&parameters); // Set defaults
  gem_mapper_parse_arguments(argc,argv,&parameters,gem_version); // Parse cmd-line

  // Runtime setup
  const mapper_parameters_cuda_t* const cuda = &parameters.cuda;
  gruntime_init(parameters.system.num_threads+1,parameters.system.tmp_folder);
  PROFILE_START(GP_MAPPER_ALL,PHIGH);
  TIMER_RESET(&parameters.loading_time);
  TIMER_RESTART(&parameters.mapper_time);

  // Open Input/Output File(s)
  TIMER_START(&parameters.loading_time);
  gem_mapper_open_input(&parameters);
  gem_mapper_open_output(&parameters);
  TIMER_STOP(&parameters.loading_time);

  // Initialize Statistics Report
  if (parameters.io.report_file_name) {
    parameters.global_mapping_stats = mm_alloc(mapping_stats_t);
  }

  // Launch mapper
  if (!cuda->gpu_enabled) {
    switch (parameters.mapper_type) {
      case mapper_se:
        mapper_se_run(&parameters);
        break;
      case mapper_pe:
        mapper_pe_run(&parameters);
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  } else {
    switch (parameters.mapper_type) {
      case mapper_se:
        mapper_cuda_se_run(&parameters);
        break;
      case mapper_pe:
        mapper_cuda_pe_run(&parameters);
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  PROFILE_STOP(GP_MAPPER_ALL,PHIGH);
  TIMER_STOP(&parameters.mapper_time);

  // Profile
  if (parameters.misc.profile) gem_mapper_print_profile(&parameters);

  // Mapping Statistics Report
  if(parameters.io.report_file_name) {
     output_mapping_stats(&parameters,parameters.global_mapping_stats);
  }

  // CleanUP
  archive_delete(parameters.archive); // Delete archive
  gem_mapper_close_input(&parameters); // Close I/O files
  gem_mapper_close_output(&parameters); // Close I/O files
  gruntime_destroy();

  // Display end banner
  const uint64_t mapper_time_sec = (uint64_t)TIMER_GET_TOTAL_S(&parameters.mapper_time);
  const uint64_t loading_time_sec = (uint64_t)TIMER_GET_TOTAL_S(&parameters.loading_time);
  gem_cond_log(parameters.misc.verbose_user,
      "[GEMMapper terminated successfully in %"PRIu64"s. (+%"PRIu64"s. loading)]\n",
      (uint64_t)BOUNDED_SUBTRACTION(mapper_time_sec,loading_time_sec,0),loading_time_sec);

  // Done!
  return 0;
}

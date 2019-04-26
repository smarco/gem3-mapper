/*
 * PROJECT: GEM-Tools library
 * FILE: gt.map2sam.c
 * DATE: 02/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Converter from MAP to SAM
 */

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

typedef struct {
  /* I/O */
  char* name_input_file;
  char* name_output_file;
  char* name_reference_file;
  char* name_gem_index_file;
  bool mmap_input;
  bool paired_end;
  bool calc_phred;
  gt_qualities_offset_t quality_format;
  /* Headers */

  /* SAM format */
  bool compact_format;
  /* Optional Fields */
  bool optional_field_NH;
  bool optional_field_NM;
  bool optional_field_XT;
  bool optional_field_XS;
  bool optional_field_md;
  /* Misc */
  uint64_t num_threads;
  bool verbose;
  /* Control flags */
  bool load_index;
  bool load_index_sequences;
} gt_stats_args;

gt_stats_args parameters = {
  /* I/O */
  .name_input_file=NULL,
  .name_output_file=NULL,
  .name_reference_file=NULL,
  .name_gem_index_file=NULL,
  .mmap_input=false,
  .paired_end=false,
  .calc_phred=false,
  .quality_format=GT_QUALS_OFFSET_33,
  /* Headers */
  /* SAM format */
  .compact_format=false,
  /* Optional Fields */
  .optional_field_NH=false,
  .optional_field_NM=false,
  .optional_field_XT=false,
  .optional_field_XS=false,
  .optional_field_md=false,
  /* Misc */
  .num_threads=1,
  .verbose=false,
  /* Control flags */
  .load_index=false,
  .load_index_sequences=false
};

gt_sequence_archive* gt_filter_open_sequence_archive(const bool load_sequences) {
  gt_sequence_archive* sequence_archive = NULL;
  if (parameters.name_gem_index_file!=NULL) { // Load GEM-IDX
    sequence_archive = gt_sequence_archive_new(GT_BED_ARCHIVE);
    gt_gemIdx_load_archive(parameters.name_gem_index_file,sequence_archive,load_sequences);
  } else {
    gt_input_file* const reference_file = gt_input_file_open(parameters.name_reference_file,false);
    sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      gt_fatal_error_msg("Error parsing reference file '%s'\n",parameters.name_reference_file);
    }
    gt_input_file_close(reference_file);
  }
  return sequence_archive;
}

#define PHRED_KONST -0.23025850929940456840 // -log(10)/10;

void gt_map2sam_calc_phred(gt_template *template)
{
	typedef struct {
		char *key;
		double prob;
		uint64_t score;
		uint8_t phred;
		UT_hash_handle hh;
	} map_hash;

	map_hash *mhash[3]={0,0,0}, *mp_hash, *tmp;
	int rd;
	size_t buf_len=1024;
	char *buf=malloc(buf_len);
	gt_cond_fatal_error(!buf,MEM_HANDLER);
	uint64_t min_score[3]={0xffff,0xffff,0xffff};
	{
		GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,maps,maps_attr) {
			if(!maps_attr) gt_fatal_error(TEMPLATE_NOT_SCORED);
			uint64_t score=maps_attr->gt_score;
			if(score==GT_MAP_NO_GT_SCORE) gt_fatal_error(TEMPLATE_NOT_SCORED);
			uint64_t seq_like[2],interval_like;
			seq_like[0]=score&0xffff;
			seq_like[1]=(score>>16)&0xffff;
			interval_like=(score>>32)&0xff;
			// Build up list of single end alignments (need this so we can scale the MAPQ score)
			// Use hash to avoid counting a single end alignment twice if it occurs in two paired alignments
			for(rd=0;rd<2;rd++) if(maps[rd]) {
				size_t ssize=gt_string_get_length(maps[rd]->seq_name);
				size_t key_size=ssize+sizeof(maps[rd]->position);
				if(key_size>buf_len) {
					buf_len=key_size*2;
					buf=realloc(buf,buf_len);
					gt_cond_fatal_error(!buf,MEM_HANDLER);
				}
				memcpy(buf,gt_string_get_string(maps[rd]->seq_name),ssize);
				memcpy(buf+ssize,&maps[rd]->position,sizeof(maps[rd]->position));
				HASH_FIND(hh,mhash[rd],buf,key_size,mp_hash);
				if(!mp_hash) {
					mp_hash=malloc(sizeof(map_hash));
					gt_cond_fatal_error(!mp_hash,MEM_HANDLER);
					mp_hash->key=malloc(key_size);
					gt_cond_fatal_error(!mp_hash->key,MEM_HANDLER);
					memcpy(mp_hash->key,buf,key_size);
					mp_hash->score=seq_like[rd];
					HASH_ADD_KEYPTR(hh,mhash[rd],mp_hash->key,key_size,mp_hash);
					if(seq_like[rd]<min_score[rd]) min_score[rd]=seq_like[rd];
				}
			}
			if(maps[0] && maps[1]) { // True paired alignments.  Shouldn't need to check for duplicates, but we will anyway
				// seq_name should be the same for the two ends in a paired alignment, but we're not taking any chances
				size_t ssize1=gt_string_get_length(maps[0]->seq_name);
				size_t ssize2=gt_string_get_length(maps[1]->seq_name);
				size_t key_size=ssize1+ssize2+2*sizeof(maps[0]->position);
				if(key_size>buf_len) {
					buf_len=key_size*2;
					buf=realloc(buf,buf_len);
					gt_cond_fatal_error(!buf,MEM_HANDLER);
				}
				memcpy(buf,gt_string_get_string(maps[0]->seq_name),ssize1);
				memcpy(buf+ssize1,gt_string_get_string(maps[1]->seq_name),ssize2);
				memcpy(buf+ssize1+ssize2,&maps[0]->position,sizeof(maps[0]->position));
				memcpy(buf+ssize1+ssize2+sizeof(maps[0]->position),&maps[1]->position,sizeof(maps[0]->position));
				HASH_FIND(hh,mhash[2],buf,key_size,mp_hash);
				if(!mp_hash) {
					mp_hash=malloc(sizeof(map_hash));
					gt_cond_fatal_error(!mp_hash,MEM_HANDLER);
					mp_hash->key=malloc(key_size);
					gt_cond_fatal_error(!mp_hash->key,MEM_HANDLER);
					memcpy(mp_hash->key,buf,key_size);
					uint64_t sc=seq_like[0]+seq_like[1]+interval_like;
					mp_hash->score=sc;
					HASH_ADD_KEYPTR(hh,mhash[2],mp_hash->key,key_size,mp_hash);
					if(sc<min_score[2]) min_score[2]=sc;
				}
			}
		}
	}
	// Now we can calculate the single and paired end MAPQ values
	for(rd=0;rd<3;rd++) if(mhash[rd]) {
		double z=0.0;
		for(mp_hash=mhash[rd];mp_hash;mp_hash=mp_hash->hh.next) {
			mp_hash->prob=exp(PHRED_KONST*(double)(mp_hash->score-min_score[rd]));
			z+=mp_hash->prob;
		}
		for(mp_hash=mhash[rd];mp_hash;mp_hash=mp_hash->hh.next) {
			mp_hash->prob/=z;
			if(1.0-mp_hash->prob<1.0e-255) mp_hash->phred=254;
			else {
				int tp=(int)(0.5+log(1.0-mp_hash->prob)/PHRED_KONST);
				if(tp>254) tp=254;
				mp_hash->phred=tp;
			}
		}
	}
	// And now we have to enter the MAPQ values in the map structures
	{
		GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,maps,maps_attr) {
			for(rd=0;rd<2;rd++) if(maps[rd]) {
				size_t ssize=gt_string_get_length(maps[rd]->seq_name);
				size_t key_size=ssize+sizeof(maps[rd]->position);
				memcpy(buf,gt_string_get_string(maps[rd]->seq_name),ssize);
				memcpy(buf+ssize,&maps[rd]->position,sizeof(maps[rd]->position));
				HASH_FIND(hh,mhash[rd],buf,key_size,mp_hash);
				assert(mp_hash);
				maps[rd]->phred_score=mp_hash->phred;
			}
			if(maps[0] && maps[1]) { // True paired alignments.  Shouldn't need to check for duplicates, but we will anyway
				// seq_name should be the same for the two ends in a paired alignment, but we're not taking any chances
				size_t ssize1=gt_string_get_length(maps[0]->seq_name);
				size_t ssize2=gt_string_get_length(maps[1]->seq_name);
				size_t key_size=ssize1+ssize2+2*sizeof(maps[0]->position);
				memcpy(buf,gt_string_get_string(maps[0]->seq_name),ssize1);
				memcpy(buf+ssize1,gt_string_get_string(maps[1]->seq_name),ssize2);
				memcpy(buf+ssize1+ssize2,&maps[0]->position,sizeof(maps[0]->position));
				memcpy(buf+ssize1+ssize2+sizeof(maps[0]->position),&maps[1]->position,sizeof(maps[0]->position));
				HASH_FIND(hh,mhash[2],buf,key_size,mp_hash);
				assert(mp_hash);
				maps_attr->phred_score=mp_hash->phred;
			}
		}
	}
	free(buf);
	for(rd=0;rd<3;rd++) if(mhash[rd]) {
		HASH_ITER(hh,mhash[rd],mp_hash,tmp) {
			HASH_DEL(mhash[rd],mp_hash);
			free(mp_hash->key);
			free(mp_hash);
		}
	}
}

void gt_map2sam_read__write() {
  // Open file IN/OUT
  gt_input_file* const input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);
  gt_output_file* const output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);
  gt_sam_headers* const sam_headers = gt_sam_header_new(); // SAM headers

  // Open reference file
  gt_sequence_archive* sequence_archive = NULL;
  if (parameters.load_index) {
    sequence_archive = gt_filter_open_sequence_archive(parameters.load_index_sequences);
    gt_sam_header_set_sequence_archive(sam_headers,sequence_archive);
  }

  // Print SAM headers
  gt_output_sam_ofprint_headers_sh(output_file,sam_headers);

  // Parallel reading+process
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    gt_status error_code;
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);

    // I/O attributes
    gt_map_parser_attributes* const input_map_attributes = gt_input_map_parser_attributes_new(parameters.paired_end);
    gt_output_sam_attributes* const output_sam_attributes = gt_output_sam_attributes_new();
    // Set out attributes
    gt_output_sam_attributes_set_compact_format(output_sam_attributes,parameters.compact_format);
    gt_output_sam_attributes_set_qualities_offset(output_sam_attributes,parameters.quality_format);
    if (parameters.optional_field_NH) gt_sam_attributes_add_tag_NH(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_NM) gt_sam_attributes_add_tag_NM(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_XT) gt_sam_attributes_add_tag_XT(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_md) gt_sam_attributes_add_tag_md(output_sam_attributes->sam_attributes);
    if (parameters.optional_field_XS) gt_sam_attributes_add_tag_XS(output_sam_attributes->sam_attributes);
    if (parameters.calc_phred) {
    	gt_sam_attributes_add_tag_MQ(output_sam_attributes->sam_attributes);
    	gt_sam_attributes_add_tag_UQ(output_sam_attributes->sam_attributes);
    	gt_sam_attributes_add_tag_PQ(output_sam_attributes->sam_attributes);
    }
    gt_template* template = gt_template_new();
    while ((error_code=gt_input_map_parser_get_template(buffered_input,template,input_map_attributes))) {
      if (error_code!=GT_INPUT_STATUS_OK) {
        gt_error_msg("Fatal error parsing file '%s':%"PRIu64"\n",parameters.name_input_file,buffered_input->current_line_num-1);
        continue;
      }
      if(parameters.calc_phred) gt_map2sam_calc_phred(template);
      // Print SAM template
      gt_output_sam_bofprint_template(buffered_output,template,output_sam_attributes);
    }

    // Clean
    gt_template_delete(template);
    gt_input_map_parser_attributes_delete(input_map_attributes);
    gt_output_sam_attributes_delete(output_sam_attributes);
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
  }

  // Release archive & Clean
  if (sequence_archive) gt_sequence_archive_delete(sequence_archive);
  gt_sam_header_delete(sam_headers);
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}

void usage(const gt_option* const options,char* groups[],const bool print_inactive) {
  fprintf(stderr, "USE: ./gt.map2sam [ARGS]...\n");
  gt_options_fprint_menu(stderr,options,groups,false,print_inactive);
  /*
   * Pending ...
   *        --score-alignments|s"
   */
  //                  "           --RG \n" // TODO: Bufff RG:Z:0 NH:i:16 XT:A:U
  // "           --headers [FILE] (Only {@RG,@PG,@CO} lines)\n"
  // "           --sorting 'unknown'|'unsorted'|'queryname'|'coordinate'\n"
}

void parse_arguments(int argc,char** argv) {
  struct option* gt_map2sam_getopt = gt_options_adaptor_getopt(gt_map2sam_options);
  gt_string* const gt_map2sam_short_getopt = gt_options_adaptor_getopt_short(gt_map2sam_options);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_map2sam_short_getopt),gt_map2sam_getopt,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 'r':
      parameters.name_reference_file = optarg;
      parameters.load_index = true;
      break;
    case 'I':
      parameters.name_gem_index_file = optarg;
      parameters.load_index = true;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 'Q':
    	parameters.calc_phred = true;
    	break;
    case 200:
      parameters.mmap_input = true;
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    /* Headers */
      // TODO
    /* Alignments */
    case 'q':
      if (gt_streq(optarg,"offset-64")) {
        parameters.quality_format=GT_QUALS_OFFSET_64;
      } else if (gt_streq(optarg,"offset-33")) {
        parameters.quality_format=GT_QUALS_OFFSET_33;
      } else {
        gt_fatal_error_msg("Quality format not recognized: '%s'",optarg);
      }
      break;
    /* Optional Fields */
    case 500: // NH
      parameters.optional_field_NH = true;
      break;
    case 501: // NM
      parameters.optional_field_NM = true;
      break;
    case 502: // XT
      parameters.optional_field_XT = true;
      break;
    case 503: // XS
      parameters.optional_field_XS = true;
      break;
    case 504: // md
      parameters.optional_field_md = true;
      break;
    /* Format */
    case 'c':
      parameters.compact_format = true;
      break;
      /* Misc */
    case 'v':
      parameters.verbose = true;
      break;
    case 't':
#ifdef HAVE_OPENMP
      parameters.num_threads = atol(optarg);
#endif
      break;
    case 'h':
      usage(gt_map2sam_options,gt_map2sam_groups,false);
      exit(1);
      break;
    case 'H':
      usage(gt_map2sam_options,gt_map2sam_groups,true);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_map2sam_options,gt_map2sam_groups,true,false);
      exit(1);
      break;
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters check
   */
  if (parameters.load_index && parameters.name_reference_file==NULL && parameters.name_gem_index_file==NULL) {
    gt_fatal_error_msg("Reference file required");
  }
  if(!parameters.load_index && parameters.optional_field_XS){
    gt_fatal_error_msg("Reference file required to compute XS field in SAM");
  }
  // Free
  gt_string_delete(gt_map2sam_short_getopt);
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  // map2sam !!
  gt_map2sam_read__write();

  return 0;
}


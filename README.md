# GEM-Mapper (Version 3)

## 1. INTRODUCTION

### 1.1. WHAT IS GEM?

GEM3 is a high-performance mapping tool for aligning sequenced reads against large reference genomes (e.g. human genome). In particular, it is designed to obtain best results when mapping sequences up to 1K bases long. GEM3 indexes the reference genome using a custom FM-Index design, and performs an adaptive gapped search based on the characteristics of the input and the user settings. GEM3 allows performing complete searches that find all possible matches according to the user criteria (i.e. avoids missing results unlike other heuristic mapping algorithms). It enables de user to perform searches looking for the first-match/best-match, all-first-n matches, and all-matches. GEM3 supports single-end, and paired-end mapping modes. Also, it supports both global-alignment and local-alignment models for different error models (i.e. hamming, edit, gap-affine). Furthermore, it offers out-of-the-box multithreaded mode to exploit several cores processors and achieve greater speed, without affecting the relative order of the reads between the input FASTQ and the ouput SAM. GEM3 is distributed under GPLv3, and runs on command line under Mac OSX and linux.

### 1.2. GETTING STARTED

```
git clone --recursive https://github.com/smarco/gem3-mapper.git gem3-mapper
cd gem3-mapper
./configure
make
```

Indexing Reference genome
```
./bin/gem-indexer -i hsapiens_v37.fa -o hsapiens_v37
```
Single-end mapping
```
./bin/gem-mapper -I hsapiens_v37.gem -i sample.Illumina.se.fastq -o sample.Illumina.se.sam
```
Paired-end mapping (single interleaved FASTA/FASTQ file)
```
./bin/gem-mapper -I hsapiens_v37.gem -i sample.Illumina.pe.fastq -p -o sample.Illumina.pe.sam
```
Paired-end mapping (splitted FASTA/FASTQ files)
```
./bin/gem-mapper -I hsapiens_v37.gem -1 sample.Illumina.pe.1.fastq -2 sample.Illumina.pe.2.fastq -o sample.Illumina.pe.sam
```

## 2. COMMAND-LINE AND OPTIONS

### 2.1 [GEM-INDEXER]
  #### I/O
```
    -i, --input=FILE (Multi-FASTA)
      Input Multi-FASTA file

    -o, --output=FILE
      Output GEM-index file prefix (.gem will be added)
```
  #### Index
```
    -b, --bisulfite-index
      Produces bisulfite compliant index
      [default=disabled]
```
  #### Performance
```
    -t, --threads=INTEGER 
      Number of threads to parallelize the mapping execution
      [default=<target-logical-cores>]
```
  #### Miscellaneous
```
    -v, --verbose={'true'|'false'}
      Enable verbose mode
      [default=true]

    --version
      Print the version number
  
    --help|-h
      Print command-line options and help
```

### 2.2 [GEM-MAPPER]

  #### I/O options
```
    -I, --index=FILE
      Input GEM index path (e.g. hsapiens_v37.gem)

    -i, --input=FILE
      Input FASTA/FASTQ file [default=stdin]
      For Single-End mode and Paired-end if the FASTA/FASTQ 
      file contains read/1 and read/2 interleaved

    -1, --i1=FILE (paired-end, end-1)
      Input Separated Paired-End FASTA/FASTQ (end/1)

    -2, --i2=FILE (paired-end, end-2)
      Input Separated Paired-End FASTA/FASTQ (end/2)

    -z, --gzip-input=FILE
      Input FASTA/FASTQ gzip compressed file

    -j, --bzip-input=FILE
      Input FASTA/FASTQ bzip compressed file

    -o, --output=FILE
      Output file [default=stdout]

    --gzip-output
      Gzip compresses the output

    --bzip-output
      Bzip compresses the output

    --report-file=FILE 
      Generates statistical report file [default=disabled]
```
  #### Single-end Alignment
```
    --mapping-mode={'fast'|'sensitive'|'customed'}
      Specify the mapping alignment approach used by the tool which ultimately 
      leads to a different accuracy/performance trade-off 
      [default=fast]

    -e, --alignment-max-error={FLOAT|INTEGER} 
      Maximum divergency rate allowed between the input sequence and the reported
      matches (i.e. maximum number of mismatches or indels allowed). All 
      global-matches above this error rate will be discarded 
      [default=0.12, 12%]

    --alignment-global-min-identity={FLOAT|INTEGER}
      Minimum global-alignment identity required (i.e. total  number of matching 
      bases). All global-matches below this identity will be discared 
      [default=80%]

    --alignment-global-min-score={FLOAT|INTEGER} 
      Minimum global-alignment score required (i.e. gap-affine score). All 
      global-matches with score below this threshold will be discarded 
      [default=0.20, 0.20*read_length*match_score]

    --alignment-local={'if-unmapped'|'never'} 
      Select whether the mapping algorithm should search for local alignments in
      case no global alignment is found, or never resort to local alignment search
      [default=if-unmapped]

    --alignment-local-min-identity={FLOAT|INTEGER} 
      Minimum local-alignment identity required (i.e. total number of matching 
      bases). All local-matches below this identity will be discarded
      [default=40]

    --alignment-local-min-score={FLOAT|INTEGER}  [default=20]
      Minimum global-alignment score required (i.e. gap-affine score). All 
      global-matches with score below this threshold will be discarded 
      [default=20]
```
  #### Paired-end Alignment
```
    -p, --paired-end-alignment
      Enable Paired-End mapping mode (in case a single FASTA/FASTQ input is given)

    -l, --min-template-length=INTEGER 
      Minimum template length allowed (distance between the beginning of the 
      left-most mapped read and the end of the right-most mapped read). All 
      matches with shorter template length will be discarded as discordant.
      [default=disabled, no restriction]

    -L, --max-template-length=INTEGER 
      Maximum template length allowed (distance between the beginning of the 
      left-most mapped read and the end of the right-most mapped read). All 
      matches with longer template length will be discarded as discordant.
      [default=10000]

    --discordant-pair-search={'always'|'if-no-concordant'|'never'} 
      Select whether the mapper should always output discordant paired matches, 
      only if no concordant paired matches are found, or never
      [default=if-no-concordant]
```
  #### Alignment Score
```
    --gap-affine-penalties=A,B,O,X 
      Gap-affine scores used by the mapping algorithm
      [default=1,4,6,1]

      A Match Score
      B Mismatch Penalty
      O Gap-Openning Penalty
      X Gap-Extension Penalty
```
  #### Reporting
```
    -M, --max-reported-matches={INTEGER|'all'} 
      Maximum number of matches reported from all the matches found by the 
      algorithm using the user settings (i.e best matches from the MAPQ-score 
      sorted list of matches found by the mapper)
      [default=5]
```
  #### Output-format
```
    -F, --output-format={'MAP'|'SAM'} 
      Select the ouput format
      [default=SAM]

    --sam-compact={'true'|'false'} 
      Output subdominant matches using the alternative hits tag in SAM format (i.e. XA)
      [default=true]

    -r, --sam-read-group-header=STRING
      Set read group header line (e.g. '@RG\tID:xx\tSM:yy')
      [default=none]
```
  #### Performance
```
    -t, --threads=INTEGER 
      Number of threads to parallelize the mapping execution
      [default=<target-logical-cores>]

    --gpu
      Enable GPU acceleration
      [default=disabled]
```
  #### Miscellaneous
```
    -v, --verbose={'quiet'|'user'|'dev'}
      Enable verbose mode 
      [default=user]

    --version
      Print the version number
  
    --help|-h
      Print command-line options and help
```
## 3. AUTHORS

  Santiago Marco-Sola \- santiagomsola@gmail.com     (Main developer)

  Alejandro Chacon    \- alejandro.chacond@gmail.com (GPU Main developer)

  Paolo Ribeca        \- paolo.ribeca@gmail.com
  
  Simon Heath         \- simon.heath@gmail.com       (Bisulfite features)

Special contributors: Jordi Camps   

## 4. REPORTING BUGS

Feedback and bug reporting it's highly appreciated. 
Please report any issue or suggestion on github, or by email to the main developer (santiagomsola@gmail.com) 

## 5. LICENSE AND CITATION

GEM3 is distributed under GPLv3 licence.

## 6. HISTORY

The GEM project was started at the CRG (Centre for Genomic Regulation, Barcelona) in 2008 by Paolo Ribeca at Guigó's Lab, where the first version of GEMv0 was developed. Later, in 2011, the project moved the development to the CNAG (Nacional Centre for Genomic Analysis, Barcelona) and Santiago Marco-Sola joined the team. During this period, GEMv1 and GEMv2 were developed, and the latter was finally released and published. In 2013, Alejandro Chacon joined the project with the purpose of developing a GPU enhanced version of the mapper (GEM-GPU). Current version (GEMv3) is actively maintained and improved by Santiago Marco-Sola at the UAB (Universitat Autonoma de Barcelona).

# gem3-mapper test suite

Tests for the mapping path of `gem-mapper`: index loading, input parsing, the
approximate search, candidate verification and alignment, match handling and
MAPQ, paired-end and bisulfite mapping, SAM/MAP output, threads, and every
command-line option. `gem-indexer` is used only to build test indexes.

There are two layers, run by one pytest command:

- **C unit tests** (`tests/c/`): small drivers that link against the objects
  of a gem3-mapper build and call internal functions directly, each checked
  against an independent reference implementation in the driver.
- **End-to-end CLI tests** (`tests/cli/`): run `gem-indexer` and `gem-mapper`
  on synthetic genomes and simulated reads whose true origin is known, and
  check the SAM output against that truth, the reference sequence and the SAM
  specification.

`TARGETS.md` lists every function `gem-mapper` reaches and every option, split
into four areas.

## Running

```bash
cd tests
python3 -m pytest                       # build gem3-mapper out of tree, run both layers
GEM3_TEST_THREADS=4 python3 -m pytest   # gem-mapper -t 4 where a test does not pin -t
python3 -m pytest c                     # C drivers only
python3 -m pytest cli -k pilot          # a subset
python3 -m pytest -k test_c_            # one area (a-d): its C drivers and CLI tests
```

Requirements: Python 3.8 or later with pytest (nothing else from PyPI), and
what gem3-mapper itself needs to build (a C compiler, make, zlib and bzip2
headers). samtools and pysam are not needed; a test that wants them must skip
without them (`pytest.importorskip("pysam")`, `shutil.which("samtools")`).
Run with `LC_ALL=C` if the locale matters to you; the harness sets it for every
process it starts.

With uv, any Python version:

```bash
cd tests && uv run --python 3.8 --with pytest python -m pytest
```

### Environment variables

| Variable | Meaning | Default |
|---|---|---|
| `GEM3_TEST_THREADS` | `-t` for tests that do not pin it | `1` |
| `GEM3_BIN_DIR` | test these binaries instead of building | build this tree |
| `GEM3_BUILD_DIR` | a configured and built gem3 tree; C drivers link its `build/*.o` (with `GEM3_BIN_DIR`) | none: C drivers skip |
| `GEM3_TEST_CACHE` | where builds and C driver binaries are cached | `$XDG_CACHE_HOME/gem3-mapper-tests` or `~/.cache/gem3-mapper-tests` |
| `GEM3_CC` | compiler given to `./configure` as `CC` | configure's choice (`gcc`) |
| `GEM3_MAKE_TARGET` | the Makefile target: `release` or `debug` (`devel` and `profile` exist but were not tried) | `release` |
| `GEM3_MAKE_ARGS` | extra make arguments, e.g. `FLAGS_OPT=-O3` | none |
| `GEM3_BUILD_JOBS` | parallel compile jobs | `4` |

### The build

Unless `GEM3_BIN_DIR` is set, the first test that needs gem3-mapper builds it
**out of tree**: the sources (everything except `.git`, `tests/` and build
products) are copied to `<cache>/build-<key>`, where the key hashes the
sources, the compiler and its version, the make target and the make
arguments, and `./configure --disable-cuda && make` runs there with the
Makefile's own flags. The source tree is never written to (a self-test checks
this). Concurrent runs wait on a lock file; a finished build is reused until a
source file changes.

The Makefiles race under `make -jN` (the tools link while objects are still
compiling, so the link fails with undefined references). The harness compiles
with `make -k -jN`, deletes whatever got linked, then links with a serial
`make`; the binaries are byte-identical to a serial build. A clean build takes
about 10 s (gcc 13, release, 4 jobs), most of it the serial LTO link.

`python3 tests/cli/gem3harness/build.py` builds (or finds) the binaries and
prints their directory. The cache is never pruned: each source state adds
about 12 MB (26 MB for a debug build), so point `GEM3_TEST_CACHE` at a
per-checkout directory when testing many versions. Because every source
change gives a new key, the first run after a change always builds, and the
self-test that the build left the source tree untouched then checks a real
build.

Indexes are built per pytest session into its temporary directory. On this
machine `gem-indexer` takes 0.3 s for the default 240 kbp genome (0.5 s for a
bisulfite index) and 2.4 s for a 2.6 Mbp genome at `-t 8` (4.9 s bisulfite),
so tests use genomes of a few hundred kbp and build extra genomes only when a
test needs a different reference.

## Layout

```
tests/
  README.md, TARGETS.md, pytest.ini, conftest.py (session fixtures)
  c/
    gem3_test.h            assert macros, TAP runner, SplitMix64 generator
    test_<name>.c          one driver per file
    test_drivers.py        builds and runs every test_*.c (item test_driver[test_<name>])
  cli/
    gem3harness/           helpers (standard library only)
      build.py             out-of-tree build, C driver build
      runner.py            Gem3: build_index(), index(), run_map(); run_process()
      genome.py            synthetic_genome(), Genome, revcomp()
      simulate.py          simulate_se(), simulate_pe(), simulate_bisulfite(), Truth
      fastx.py             FASTQ/FASTA writers: plain, gzip, bzip2; interleaved, split
      sam.py               SAM parser, CIGAR/MD/NM helpers, sam_digest()
      rng.py               Rng (SplitMix64), the same stream on every platform
    test_harness_selftest.py
    test_pilot_mapping.py
```

## Adding a test

**A CLI test** is a function in a `tests/cli/test_*.py` file using the session
fixtures from `conftest.py`:

| Fixture | What it is |
|---|---|
| `gem3` | `gem3harness.Gem3` for the build under test |
| `gem3_build` | the build (`bin_dir`, `root`, `info` with compiler and flags) |
| `genome` | `synthetic_genome(seed=1)`: ~240 kbp, chr1-3 plus 400 bp and 90 bp contigs |
| `genome_index`, `genome_bs_index` | its regular and bisulfite (`-b`) `.gem` paths |
| `threads` | `GEM3_TEST_THREADS` |
| `tmp_path` | pytest's per-test directory |

```python
from gem3harness import simulate_se, write_reads, check_record_against_reference

def test_example(gem3, genome, genome_index, tmp_path, threads):
	reads = simulate_se(genome, 200, 100, seed=3, subs=1)
	fq = write_reads(tmp_path / "r.fq", reads)
	res = gem3.run_map(genome_index, inputs=[fq], args=["-M", "1"], threads=threads)
	sam = res.sam()
	for rec in sam.records:
		assert check_record_against_reference(rec, genome[rec.rname]) == []
```

Main helpers (docstrings have the details):

- `synthetic_genome(seed, contigs=..., ...)` -> `Genome` (`genome[name]` is the
  sequence with case kept; `genome.features` lists repeats, tandems,
  homopolymers, N runs, soft-masked stretches and short contigs;
  `genome.hard_regions(contig, start, end)`; `write_fasta(path)`).
- `simulate_se(genome, n, length, seed, subs/ins/dels or *_rate, unique_only, strand, ...)`,
  `simulate_pe(..., insert_mean, insert_sd, orientation="FR"|"RF"|"FF"|"RR", adapters=...)`,
  `simulate_bisulfite(..., library="directional"|"pbat"|"non-directional", cpg_rate, chh_rate, paired)`.
  Each read has `.truth`: `contig`, `pos` (0-based), `strand`, `cigar`, `md`,
  `nm`, `ref_span`, `genomic_length`, `edits`, `extra` (adapter length,
  bisulfite origin OT/OB/CTOT/CTOB, conversion, methylation calls).
- `write_reads(path, reads, fmt="fastq"|"fasta", compression=None|"gz"|"bz2", suffix="")`,
  `write_interleaved(path, pairs, ...)`, `write_split(p1, p2, pairs, ...)`;
  mate names get `/1` and `/2`.
- `gem3.index(genome, bisulfite=False)` / `gem3.build_index(fasta_path, ...)`:
  cached per FASTA bytes and options for the session.
- `gem3.run_map(index, args=(), inputs=(), stdin=None, threads=None, timeout=120, check=True, quiet=True, files=None)`
  -> `RunResult` (`returncode`, `stdout`, `stderr` bytes, `stderr_text`,
  `files` written to the working directory, `sam()`, `timed_out`,
  `crashed`). It adds `--verbose=quiet` unless `quiet=False` or `-v` is in
  `args`, runs in a private temporary directory, with `LC_ALL=C`, core dumps
  off and a timeout. `check=False` for tests that expect a failure.
- `parse_sam(text)` -> `Sam` (`header_lines`, `sq`, `header("@RG")`,
  `records`, `by_name()`, `primary()`); `SamRecord` fields with typed `tags`.
- `md_nm_from_alignment(ref, pos, cigar_ops, seq)`, `bisulfite_nm(...)`,
  `cigar_ref_span`, `cigar_query_length`, `cigar_is_valid`,
  `check_record_against_reference(rec, ref_seq)`.
- `sam_digest(data)` / `strip_variable_header(data)`: compare or pin SAM bytes
  without the `@PG` and `@CO` lines.

**A C driver** is one new file, `tests/c/test_<name>.c`; pytest finds it. It
starts with gem3's GPL header block, includes `gem3_test.h` (which includes
gem3's `utils/essentials.h` first) and the gem3 headers it needs, and ends with
the case list, which defines `main`:

```c
#include "gem3_test.h"
#include "align/align_bpm_distance.h"

GT_TEST(my_case) {
  gt_rng_t rng = gt_rng_new(7);
  GT_CHECK_EQ_UINT(1+1,2);
  GT_REQUIRE(gt_rng_below(&rng,10) < 10);
}

GT_CASES(
  GT_CASE(my_case)
)
```

`GT_CHECK`, `GT_CHECK_MSG(cond, fmt, ...)`, `GT_CHECK_EQ_INT/UINT/STR/MEM`
record a failure and continue; `GT_REQUIRE*` end the case. A driver that
starts threads runs them with `gt_run_threads(n, fn, args, arg_size)`, which
gives them gem thread ids 1..n and the profile slots that `make debug`
(`GEM_PROFILE`) builds count into, as gem-mapper does for its mapper
threads; check results after it returns, not inside `fn`. A driver prints TAP
and exits non-zero on failure; run one by hand with
`<cache>/c-drivers/test_<name>-<key>` (`-l` lists cases, names select them). It
is compiled with the flags gem3's Makefile used for `gem-mapper` (for the
default build `-Ofast -march=skylake-avx512 -flto` on this kind of CPU) and
linked with every object in the build's `build/` directory. A driver that
calls an internal function says in its header comment which internal
signatures it uses and which struct fields it reads or sets, since a speed
change may change them. A change to those signatures or fields may edit a
driver at those call sites and field accesses only, never its expected
values or inputs.

Drivers that build gem3's internal structs field by field
(`test_c_matches.c`, `test_c_paired.c`, `test_d_stats.c`) need more than
call-site edits if a change reorganises those structs (merges fields,
splits a struct, changes what a field means); such a change should keep the
driver's checks and move them to the new layout, or replace them with a CLI
test of the same behaviour.

## Rules for expected values

These bind every test in the suite.

- Derive expected values independently wherever practical: inputs small
  enough to work by hand; a short reference implementation inside the test
  (brute-force edit distance, gap-affine alignment with gem3-mapper's scoring
  `--gap-affine-penalties` A,B,O,X = match, mismatch, gap open, gap extension,
  where a gap of length k costs O + k*X; NM and MD recomputed from the genome;
  a naive FM-index count on a tiny text); the simulator's recorded truth; the
  SAM specification. A test that only records the current output says so in
  its docstring ("pins the current output"); keep these few.
- No test of behaviour the base version gets wrong. Such behaviour is a known
  bug, recorded outside the suite with a reproducer (see "Known bugs" below).
- No dependence on platform, compiler, build flags, Python version or thread
  count, unless the test is explicitly about thread-count invariance.
- The whole suite (both layers, build excluded) stays under about 2 minutes
  at 4 threads on a loaded machine; it runs once per change during speed work.
- Do not pin implementation details that a speed optimization could change
  while keeping `gem-mapper`'s output byte-identical: internal struct layouts,
  allocation counts, the order of internal calls, buffer sizes, private
  function names. Pinning the bytes of `gem-mapper` SAM output on fixed inputs
  is fine and wanted (use `sam_digest`). C tests of internal functions test
  behaviour through inputs and outputs and name the internal signatures they
  call.
- Fixtures are generated by the tests (seeded, `gem3harness.rng.Rng`, never
  `random`'s `choice`/`randrange`/`gauss`, whose streams may change between
  Python versions). Small subsets of public reads (ENCODE ENCSR747ZAI /
  ENCFF459LYK ChIP-seq, SRA SRR13051253 WGBS) may be added with the accession
  in a README beside them; GIAB HG002 WGS reads may not.

## Determinism (measured on base 6e155e3, 2026-10-29)

Inputs: a 2.6 Mbp synthetic genome (seed 7, repeats included) with 60,000
single-end 100 bp reads, 50,000 2x100 pairs (3% with a mate from elsewhere)
and 20,000 directional bisulfite pairs; plus 38 option configurations on 1,000-read
subsets (mapping modes, -M, MAPQ model, MAP output, `--sam-compact=false`,
alignment models and penalties, local alignment, candidate generation and
verification, template-length limits, discordant search, bisulfite read
modes).

- **Threads and repeats.** SAM records and every header line except `@PG` and
  `@CO` are byte-identical at `-t 1, 4, 8, 32` and across repeated runs, SE,
  PE and bisulfite. Output order follows the input order at every `-t`. This
  holds for the default `--input-model` on inputs smaller than (processors x
  32 MB), where no input block is refilled. Blocks that get refilled race at
  `-t` above 1 (reads lost or garbled in about 2% of runs with 512 B-2 KB
  blocks, or a hang; `known-bugs/r1-input-model-reader-count-race.md`), and
  blocks of a few KB drop lines
  (`known-bugs/A-input-model-small-blocks-drop-lines.md`); the default model
  refills only for inputs above 8 GB on this 256-core machine.
- **Header lines that vary.** `@PG` holds `VN:` (the version string: the
  `VERSION` file plus `-release`/`-debug`, or `git describe` when the build tree
  has a `.git`) and `CL:` (the command line after getopt permuted argv and after
  comma-separated option values were split in place, so `--input-model
  1M,2,64K` shows as `--input-model 1M`). `@CO` holds `TM:` (local time to the
  second), `WD:` (working directory), `HN:` (host) and `UN:` (user). Compare
  with `sam_digest`, which drops both lines.
- **Per-buffer state.** Paired-end output depends on the input buffer size
  (`--input-model`'s third field, default 4M): the template-length estimate is
  reset at every input-buffer reload (`mapper_pe_reload_buffer` in
  `src/mapper/mapper_io.c`), so pairs whose placement is a tie between
  equal-scoring candidates (tandem repeats, homopolymers) can resolve
  differently with other buffer sizes. Because each buffer is mapped whole by
  one thread, the thread count does not matter. Single-end output did not
  change with any buffer or block size tried (256K to 64M blocks, 16K to 16M
  buffers). A speed change must keep the input buffer boundaries and this
  reset.
- **Builds.** All 38 configurations and the three larger inputs give
  byte-identical SAM (without `@PG`/`@CO`) from: gcc 13.3 release (`-Ofast
  -march=skylake-avx512 -flto`, the Makefile's choice here), gcc 13.3 `make
  debug` (`-Og -g -DGEM_DEBUG -DGEM_PROFILE`), gcc 13.3 with `FLAGS_OPT=-O3`,
  gcc 11.5 with the release flags, gcc 11.5 as `CC=gcc-11` (the Makefile then
  falls back to `-O3`: its optimization flags apply only when `CC` is exactly
  `gcc`), and clang 18.1 (`-O3`). Every build writes a byte-identical index.
  These inputs miss one case where builds disagree
  (`known-bugs/r2-mapq-infinite-to-integer.md`): a read whose MAPQ logistic
  probability rounds to exactly 1.0 (linear predictor above about 36.7, e.g.
  100 bp reads missing 12-18 bp, a 250 bp read under `--alignment-model
  hamming`, some reads under `--candidate-generation MEM`) has
  `round(-10*log10(0))` = +inf converted to an integer, which C leaves
  undefined. It prints MAPQ 59 (classify 170 or 70) from AVX-512 gcc builds
  and clang, and 0 from gcc builds without AVX-512 (`-O3`, `-Og`,
  `CC=gcc-11`, `-Ofast -march=native` on a CPU without AVX-512). It is not
  floating-point contraction or fast-math: the predictors printed by
  `--mapq-model dump-predictors` and the probability itself are identical in
  every build. No test checks MAPQ on a read whose probability rounds to
  1.0 (area A's MEM pin in `test_a_candidate_generation.py` leaves MAPQ out,
  and area B's pins leave those reads out); do not pin one. The search's
  reference is
  the Makefile's own build (gcc 13 `-Ofast -march=skylake-avx512 -flto` on
  this machine), with the same flags for the baseline and every candidate.
- **Uninitialized and stale memory** (checked with valgrind memcheck on a
  `-O0 -g` build, `MALLOC_PERTURB_`, ASan, and a build that logs the reads in
  question; 45 option configurations). Output depends on memory gem-mapper
  never wrote in three cases, which the search's gate cases must avoid:
  `--mapping-mode complete-partition` and `complete-partition-filter`, whose
  MAPQ comes from match metrics that are never set (read length 0 in a
  zeroed heap gives NaN predictors, which every x86 build here converts to
  59; ASan's heap fill turns 26 and 128 of area B's records to MAPQ 0;
  `known-bugs/r1-complete-partition-mapq-uninitialized.md`; area B's pins of
  these modes leave MAPQ out);
  paired-end `--sam-compact=false` at any `-M` (`XS:i` indexes one end's
  alignments by pair number, so a template with more pairs than that end has
  alignments prints a stale alignment from an earlier read of the same thread
  or crashes; `known-bugs/D-pe-sam-compact-false-crash.md`), and paired-end
  bisulfite mapping whatever `--bisulfite-read` says (the paired-end I/O
  handler never sets `bisulfite_read`; directional pairs gave the same bytes
  with the field forced to `non-stranded`, PBAT and non-directional pairs did
  not; `known-bugs/pe-bisulfite-non-stranded-ignored.md`); the suite's
  paired bisulfite tests use directional pairs and pass with the field
  forced to `non-stranded` (an `LD_PRELOAD` malloc fill). The suite's two
  paired `--sam-compact=false` cases (`test_pair_tie_between_two_copies`, pin
  `rg-Mall-full`) have as many alignments per end as pairs; a probe build
  logged no read past an end's alignments. Indexes depend on memory too:
  gem-indexer ORs the sampled suffix array into slab memory it never zeroes,
  which is right only because fresh pages are zero
  (`known-bugs/r2-indexer-sampled-sa-uninitialized.md`). Sanitizer runs must
  build their indexes with an ordinary gem-indexer, and a change that makes
  `mm_slab` or `segmented_vector` memory come back non-zero corrupts every
  index the suite builds.
- **Tied pairs.** Among concordant pairs that tie on relation, SWG score and
  whole sigmas of template length, the primary pair (and the order of the
  secondary records and XA) is the one enumerated first:
  `paired_matches_cmp_swg_score` returns 0 between pairs whose layouts differ
  (its layout line compares the relations, equal by then), and pairs that
  share their leftmost mate tie on `index_position`
  (`known-bugs/r2-pe-tied-pair-order.md`). Reversing the enumeration moved the
  primary placement of 9 of 350 tie-prone templates. A speed change must keep
  the enumeration order of `paired_matches_find_pairs_sorted`, the order of
  the `qsort` by genomic position before it, and insertion after equal keys.
  `test_paired_output_follows_input_buffers_not_threads` and
  `test_tandem_array_pairs_are_pinned` (the latter also at `-M all`, which
  catches a change to the layout tie-break) pin such templates.
- **Sanitizers** (gcc 13 `-O1 -fsanitize=address,undefined`, areas C and D
  and the pilot tests, 2026-10-30): no AddressSanitizer error in gem-mapper
  or the C drivers. UBSan reports misaligned accesses to index structures
  used in place (the locator, `locator.c`, `mm.c:522`, `output_sam.c:105-115`;
  their offset in the file is not a multiple of 8), `1L << 63`
  (`align_bpm.c:208`) and, in gem-indexer, a shift by -2 when printing the
  text sampling rate of an index built without text sampling
  (`sampled_sa.c:153`, only the `.info` file); none changes output on
  x86-64 (`known-bugs/r2-sanitizer-findings.md`). A sanitizer run of the
  suite needs `-fno-sanitize=alignment`, `ASAN_OPTIONS=max_malloc_fill_size=0`
  (so gem-indexer gets zeroed memory) and `detect_leaks=0`.
- **Possibly platform-dependent, not testable here**: MAPQ under `-ffast-math`
  on aarch64/ppc64le (`-Ofast -mcpu=native`); the SIMD SWG paths
  (`align_swg_simd.c`) selected by `-march`; `long`/`size_t` widths on other
  ABIs; macOS's getopt (BSD getopt does not permute argv, so the `@PG CL`
  order and the handling of stray arguments may differ); `-lrt` and
  `/proc/cpuinfo`-based flag selection in `Makefile.mk.in`.

### Compiler setups used here

```bash
cd tests
python3 -m pytest                                          # gcc 13, release flags
GEM3_MAKE_TARGET=debug python3 -m pytest                   # -Og, GEM_DEBUG checks
GEM3_MAKE_ARGS=FLAGS_OPT=-O3 python3 -m pytest             # portable flags
GEM3_CC=gcc-11 python3 -m pytest                           # gcc 11 (Makefile uses -O3)
PATH=<dir with gcc -> gcc-11, gcc-ar -> gcc-ar-11>:$PATH python3 -m pytest   # gcc 11, release flags
GEM3_CC=<clang wrapper> python3 -m pytest                  # clang 18 (-O3)
```

A conda-forge clang needs the host's headers and libraries for zlib and
bzip2; a wrapper script calling it with `--target=x86_64-linux-gnu
--sysroot=/ --gcc-toolchain=/usr "$@"` works.

## gem3-mapper behaviour test writers should know

- `-v`, `--sam-compact`, `--clipping`, `--alignment-force-full-swg`,
  `--alignment-curation`, `--sam-gem-compatible`, `--profile` and `-h` take an
  optional value, which must be attached: `--sam-compact=false`,
  `--verbose=quiet`. `--sam-compact false` silently means `--sam-compact`.
  Stray positional arguments are ignored without a message.
- Default output: `@HD VN:1.4 SO:unsorted`, one `@SQ` per contig in FASTA
  order (a bisulfite index adds no contigs), `@PG`, `@CO`. Mapped records carry
  `NM:i`, `MD:Z` (not on bisulfite reads), `AS:i` (the gap-affine score of the
  reported alignment; default penalties 1,4,6,1, a k-base gap costs 6 + k) and,
  for multi-mapping reads in compact mode, `XA:Z`
  (`contig,±pos,CIGAR,NM;`). With `--sam-compact=false` the alternatives are
  separate records with flag 256.
- SEQ and QUAL are reverse-complemented/reversed for reverse-strand records.
  Read names lose a `/1` or `/2` suffix.
- Bisulfite (`-b` index): `XB:A:C` when the read aligns to the C-to-T
  converted forward strand, `XB:A:G` for G-to-A; `NM` does not count
  conversions (`bisulfite_nm` recomputes it); SEQ is the read as sequenced.
  `--bisulfite-read inferred` (the default) takes the conversion from the read
  name's `/1` or `/2` suffix in single-end mode; reads without a suffix are not
  converted. In paired-end mode end 1 is C-to-T and end 2 G-to-A.
- Paired-end mode is set by `-1/-2` or by `-p` (interleaved input, 8 lines per
  pair in FASTQ). TLEN is positive on the leftmost mate.
- MAPQ values seen with the default model ranged from 0 to 60; unique reads
  with several errors can get a low MAPQ (e.g. 6) while placed correctly.
- `--mapping-mode complete-brute-force` took over 60 s on 1,000 reads against
  a 2.6 Mbp genome; keep inputs tiny for it.
- `gem-indexer` writes `<name>.sa.tmp` into its working directory, not beside
  the output prefix; `Gem3.build_index` runs it in the index's own directory.

## Coverage (measured 2026-10-30)

gem3-mapper had no tests of its own (0 tests, 0% line coverage). The whole
suite (C drivers and CLI tests, `GEM3_TEST_THREADS=4`) was run against a
gcc 13 build with `FLAGS_OPT="-O0 --coverage -fprofile-abs-path
-fprofile-update=atomic"` and read with gcovr 8.8 (`gcovr --json`); the
gem-indexer runs that build the session's indexes count too. Files under
`resources/` and `src/gpu/` are left out.

- **All of `src/` and `include/`:** 13,417 of 20,861 executable lines
  (64.3%) in 191 files.
- **The mapping path** (the 138 files whose functions `TARGETS.md` lists):
  12,114 of 15,486 lines (78.2%). `TARGETS.md` ticks 232 of 252 area A
  functions, 103 of 334 in B (most of the rest run under the CLI tests and
  are marked `indirect`), 141 of 149 in C, 250 of 276 in D, and 54 of 69
  options.

| Part | Lines | Covered | Tests |
|---|---:|---:|---:|
| `src/align/` (+ `pattern/`, `include/utils/`) | 1,709 | 73.4% | 127 |
| `src/approximate_search/` | 660 | 46.7% | 252 |
| `src/archive/` | 505 | 91.3% | 125 |
| `src/archive/score/` | 242 | 97.1% | 122 |
| `src/archive/search/` | 886 | 54.6% | 249 |
| `src/archive/builder/`, `src/fm_index/sa_builder/` (gem-indexer only) | 769 | 89.2% | 0 |
| `src/filtering/candidates/` | 1,334 | 43.7% | 127 |
| `src/filtering/region/` | 253 | 77.5% | 127 |
| `src/filtering/region_profile/` | 909 | 77.2% | 125 |
| `src/fm_index/` (+ `bwt/`) | 1,210 | 63.0% | 125 |
| `src/io/` | 1,867 | 91.0% | 269 |
| `src/mapper/` | 904 | 48.3% | 144 |
| `src/matches/` | 935 | 81.1% | 122 |
| `src/matches/align/`, `scaffold/` | 1,153 | 88.2% | 127 |
| `src/matches/classify/` | 307 | 91.2% | 122 |
| `src/neighborhood_search/` | 883 | 72.5% | 127 |
| `src/profiler/` | 96 | 29.2% | 144 |
| `src/search_pipeline/` (CUDA mapper only) | 881 | 0.0% | 0 |
| `src/stats/` | 832 | 53.4% | 144 |
| `src/system/` | 1,483 | 61.1% | 269 |
| `src/text/` | 597 | 35.3% | 247 |
| `src/tools/` (+ `interface/`) | 1,249 | 63.7% | 144 |
| `src/utils/` | 1,197 | 44.3% | 396 |

"Tests" counts the pytest items and C driver cases in the files of the areas
that own functions in that directory (`TARGETS.md`), so a test can count in
several rows; the suite has 583 in all (467 pytest items under `cli/`,
harness self-tests included, and 116 cases in the 15 C drivers). Much of what the low rows miss is in files no
test executes: `src/mapper/mapper_cuda*.c` (413 of the 467 lines missed in
`src/mapper/`; the CUDA mapper), `src/text/cdna_text.c` and
`cdna_bitwise_text.c` (314 lines), and `src/utils/sparse_array_locator.c`,
`interval_set.c`, `sparse_bitmap.c` and `region_set.c` (392 lines).

Mapping-path functions that no test calls (not executed in this run and not
named in a test): the debug and `-c` failure printers
`alignment_print_pretty`, `archive_check_se_match_print` (and its
`_incorrect` and `_suboptimum` variants), `match_alignment_print_pretty`,
`match_alignment_region_print`, `match_scaffold_print`;
`nsearch_levenshtein_base` and `nsearch_levenshtein_base_step` (reached only
through a known bug of `--candidate-generation`);
`nsearch_hamming_preconditioned` (no input found that reaches it);
`filtering_candidates_sort_regions_by_scaffold_coverage`,
`filtering_candidates_classify_subdominant_match_edit`,
`filtering_candidates_extend_candidates`; `mapper_profile_print_mapper_pe`
(empty outside `GEM_PROFILE` builds); and `mm_allocator_free_malloc` (only
with `MM_ALLOCATOR_MALLOC` defined).

## Known bugs

Base behaviour that contradicts the help text, the source's own intent or the
SAM specification is not tested. Found while writing this suite:

- `--pair-orientation`, `--discordant-pair-orientation`, `--pair-layout`,
  `--discordant-pair-layout`: any valid value loops forever while parsing.
- `--pe-template-length min,max,samples` is rejected; two values crash.
- `--candidate-generation adaptive` / `adaptive-limited` are rejected.
- A space-separated value of an optional-argument option is ignored
  (upstream #37), as are stray positional arguments.
- Single-end `--bisulfite-read 2` and `interleaved` behave like `1`.
- Paired-end `--bisulfite-read non-stranded` is never applied (the paired-end
  I/O handler leaves `bisulfite_read` uninitialized), so PBAT pairs do not map.
- `--input-model`: blocks of a few KB drop lines, refilled blocks race at
  `-t` above 1, one block hangs once the input exceeds it, and two blocks
  hang paired-end at `-t 32`.
- `--mapping-mode complete-partition` and `complete-partition-filter` compute
  MAPQ from match metrics that are never set.
- Paired-end `--sam-compact=false` prints `XS:i` from another read's
  alignment or crashes when one end has fewer alignments than the template
  has pairs.
- MAPQ of a read whose logistic probability rounds to 1.0 depends on the
  compiler (+inf converted to an integer).
- Tied concordant pairs are ordered by enumeration order (the layout
  tie-break compares the wrong field).
- gem-indexer builds a correct index only from zeroed heap memory.
- Sanitizers (ASan, UBSan) report an overlapping `memcpy` in the FASTQ
  header parser, a read past `dna_encoded_complement_table` on run-length
  indexes, a use after free on a failing input pipe, misaligned loads from
  the index, undefined shifts and leaks at exit (`r1-sanitizer-findings.md`,
  `r2-sanitizer-findings.md`).

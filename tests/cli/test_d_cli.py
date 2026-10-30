# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""The gem-mapper command line (src/tools/gem-mapper.c,
src/tools/interface/mapper_arguments.c, src/utils/options_menu.c,
src/mapper/mapper_parameters.c) and its log output (src/system/report.c).

Known bugs left out (known-bugs/): -t, or the buffer and block counts of
--output-model and --input-model, given as 0 or a non-number map nothing
(D-zero-or-non-numeric-counts.md); unequal pair files and an odd
interleaved file print the error and then hang (D-unsynced-pairs-hang.md);
-B and --gap-affine-penalties flip the mismatch score's sign
(B-mismatch-penalty-sign.md), so the documented 1,4,6,1 is not given
explicitly; --alignment-curation-min-end-context and
--alignment-global-min-score have other defaults than --help states
(B-help-defaults.md); --map-format 1 aborts (D-map-format-1-not-implemented.md).
"""

import os
import pwd
import re

import pytest

import area_d
from gem3harness import run_process, sam_digest, strip_variable_header

# gem_mapper_options[]: (long option, short option or None, visibility).
OPTIONS = [
	("index", "I", "user"), ("input", "i", "user"), ("i1", "1", "user"), ("i2", "2", "user"),
	("gzip-input", "z", "user"), ("bzip-input", "j", "user"), ("input-model", None, "developer"),
	("output", "o", "user"), ("gzip-output", None, "user"), ("bzip-output", None, "user"),
	("output-model", None, "developer"), ("report-file", None, "user"), ("clipping", None, "advanced"),
	("quality-format", "q", "advanced"), ("mapping-mode", None, "user"),
	("complete-search-error", "E", "advanced"), ("complete-strata-after-best", "s", "advanced"),
	("alignment-max-error", "e", "user"), ("alignment-max-bandwidth", None, "advanced"),
	("alignment-global-min-identity", None, "user"), ("alignment-global-min-score", None, "user"),
	("alignment-local", None, "user"), ("alignment-local-min-identity", None, "user"),
	("alignment-local-min-score", None, "user"), ("alignment-force-full-swg", None, "advanced"),
	("alignment-scaffolding-min-matching_length", None, "developer"),
	("alignment-curation", None, "advanced"), ("alignment-curation-min-end-context", None, "developer"),
	("candidate-generation", None, "developer"), ("candidate-generation-adaptive", None, "developer"),
	("candidate-verification", None, "developer"), ("qgram-filter", None, "developer"),
	("paired-end-alignment", "p", "user"), ("min-template-length", "l", "user"),
	("max-template-length", "L", "user"), ("discordant-pair-search", None, "user"),
	("pair-orientation", None, "advanced"), ("discordant-pair-orientation", None, "advanced"),
	("pair-layout", None, "advanced"), ("discordant-pair-layout", None, "advanced"),
	("pe-template-length", None, "advanced"), ("bisulfite-read", None, "user"),
	("underconversion_sequence", None, "user"), ("overconversion_sequence", None, "user"),
	("control_sequence", None, "user"), ("alignment-model", None, "advanced"),
	("gap-affine-penalties", None, "user"), ("matching-score", "A", "advanced"),
	("mismatch-penalty", "B", "advanced"), ("gap-open-penalty", "O", "advanced"),
	("gap-extension-penalty", "X", "advanced"), ("mapq-model", None, "advanced"),
	("min-reported-strata", "m", "advanced"), ("max-reported-matches", "M", "user"),
	("output-format", "F", "user"), ("sam-compact", None, "user"),
	("sam-read-group-header", "r", "user"), ("sam-gem-compatible", None, "developer"),
	("map-format", None, "developer"), ("threads", "t", "user"), ("max-memory", None, "developer"),
	("tmp-folder", None, "developer"), ("check-alignments", "c", "developer"),
	("profile", None, "developer"), ("verbose", "v", "user"), ("progress-step", None, "developer"),
	("version", None, "user"), ("help", "h", "user"),
]
LEVELS = ("user", "advanced", "developer")
LOG_LINE = re.compile(r"\A\d{4}/\d{1,2}/\d{1,2} \d\d:\d\d:\d\d -- (.*)\Z")


@pytest.fixture(scope="module")
def se(gem3, genome, genome_index):
	return area_d.se_dataset(gem3, genome, genome_index)


@pytest.fixture(scope="module")
def pe(gem3, genome, genome_index):
	return area_d.pe_dataset(gem3, genome, genome_index)


def mapper(gem3, *args, **kw):
	return run_process([gem3.mapper] + [str(a) for a in args], **kw)


def help_options(text):
	"""(long, short) of every option line in a usage text, in order."""
	out = []
	for line in text.split("\n"):
		m = re.match(r"\A      --([A-Za-z0-9_-]+)(\|-(\S))? ", line)
		if m:
			out.append((m.group(1), m.group(3)))
	return out


# ---------------------------------------------------------------------------
# --help, --version, no arguments
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("level", LEVELS)
def test_help_lists_options_by_visibility(gem3, level):
	"""--help=<level> prints on stderr, in table order, the options visible
	at that level (user < advanced < developer), each under its group."""
	res = mapper(gem3, "--help=%s" % level)
	assert res.returncode == 0 and res.stdout == b""
	text = res.stderr_text
	assert text.startswith("USAGE: ./gem-mapper [ARGS]...\n")
	allowed = LEVELS[:LEVELS.index(level) + 1]
	assert help_options(text) == [(l, s) for l, s, v in OPTIONS if v in allowed]
	groups = re.findall(r"^    \[(.+)\]$", text, re.M)
	assert groups[0] == "I/O" and groups[-1] == "Miscellaneous"
	assert "CUDA Settings" not in groups


def test_help_spellings_agree(gem3):
	"""No arguments, -h, --help, --help=user and --help=USER print the same
	user-level usage and exit 0."""
	ref = mapper(gem3, "--help=user")
	for args in ([], ["-h"], ["--help"], ["--help=USER"], ["-huser"]):
		res = mapper(gem3, *args)
		assert (res.returncode, res.stdout, res.stderr) == (0, b"", ref.stderr), args


def test_help_documents_defaults(gem3):
	"""Spot checks of the default values printed by --help=developer."""
	text = mapper(gem3, "--help=developer").stderr_text
	for frag in ("--max-reported-matches|-M <number>|'all' (default=5)",
			"--output-format|-F 'MAP'|'SAM' (default=SAM)",
			"--output-model <buffer_size,num_buffers> (default=4M,5c)",
			"--input-model <block_size,num_blocks,buffer_size> (default=32M,c,4M)",
			"--gap-affine-penalties A,B,O,X (default=1,4,6,1)"):
		assert frag in text, frag


def test_version(gem3, se):
	"""--version prints the version on stderr and exits; @PG VN holds it."""
	res = mapper(gem3, "--version")
	assert res.returncode == 0 and res.stdout == b""
	v = res.stderr_text.strip()
	assert re.match(r"\Av\d+\.\d+\.\d+", v), v
	assert se.run().sam().header("@PG")[0]["VN"] == v


# ---------------------------------------------------------------------------
# Argument errors
# ---------------------------------------------------------------------------

ERRORS = [
	("no-index", ["-i", "{se}"], "Index file required"),
	("index-is-input", ["-I", "{idx}", "-i", "{idx}"], "Index and Input-file must be different"),
	("input-is-output", ["-I", "{idx}", "-i", "{se}", "-o", "{se}"], "Input-file and Output-file must be different"),
	("index-is-output", ["-I", "{idx}", "-i", "{se}", "-o", "{idx}"], "Index and Output-file must be different"),
	("missing-i2", ["-I", "{idx}", "-1", "{p1}"], "Missing Input-End2 (--i2)"),
	("missing-i1", ["-I", "{idx}", "-2", "{p2}"], "Missing Input-End1 (--i1)"),
	("same-ends", ["-I", "{idx}", "-1", "{p1}", "-2", "{p1}"], "Input-End1 and Input-End2 must be different"),
	("index-is-end1", ["-I", "{idx}", "-1", "{idx}", "-2", "{p2}"], "Index and Input-End1 must be different"),
	("index-is-end2", ["-I", "{idx}", "-1", "{p1}", "-2", "{idx}"], "Index and Input-End2 must be different"),
	("end1-is-output", ["-I", "{idx}", "-1", "{p1}", "-2", "{p2}", "-o", "{p1}"], "Input-End1 and Output-file must be different"),
	("end2-is-output", ["-I", "{idx}", "-1", "{p1}", "-2", "{p2}", "-o", "{p2}"], "Input-End2 and Output-file must be different"),
	("m-and-M", ["-I", "{idx}", "-i", "{se}", "-m", "1", "-M", "2"], "cannot be used at the same time"),
	("M-zero", ["-I", "{idx}", "-i", "{se}", "-M", "0"], "must be greater than zero"),
	("unknown", ["-I", "{idx}", "-i", "{se}", "--bogus"], "Option not recognized"),
	("gpu", ["-I", "{idx}", "-i", "{se}", "--gpu"], "Option not recognized"),
	("gpu-devices", ["-I", "{idx}", "-i", "{se}", "--gpu-devices", "all"], "Option not recognized"),
	("gpu-buffers", ["-I", "{idx}", "-i", "{se}", "--gpu-buffers-model", "2,3,3,3,3,1M"], "Option not recognized"),
	("candidate-drop-off", ["-I", "{idx}", "-i", "{se}", "--candidate-drop-off", "1"], "Option not recognized"),
	("pe-extension", ["-I", "{idx}", "-i", "{se}", "--pe-extension", "1"], "Option not recognized"),
	("output-format", ["-I", "{idx}", "-i", "{se}", "-F", "BAM"], "must be 'MAP' or 'SAM'"),
	("map-format", ["-I", "{idx}", "-i", "{se}", "--map-format", "4"], "must be in {'1'|'2'|'3'}"),
	("verbose", ["-I", "{idx}", "-i", "{se}", "--verbose=loud"], "must be 'quiet'|'user'|'dev'"),
	("profile", ["-I", "{idx}", "-i", "{se}", "--profile=median"], "must be 'sum'|'min'|'max'|'mean'|'sample'"),
	("help", ["--help=expert"], "Help argument not valid"),
	("max-memory", ["-I", "{idx}", "-i", "{se}", "--max-memory", "lots"], "not a valid size"),
	("output-model-count", ["-I", "{idx}", "-i", "{se}", "--output-model", "4M"], "'--output-model' wrong number of arguments"),
	("output-model-size", ["-I", "{idx}", "-i", "{se}", "--output-model", "x,2"], "Error parsing 'buffer_size'"),
	("rg-prefix", ["-I", "{idx}", "-i", "{se}", "-r", "ID:x"], "must start with '@RG'"),
	("rg-id", ["-I", "{idx}", "-i", "{se}", "-r", "@RG\\tSM:x"], "must contain ID field"),
	("rg-long-id", ["-I", "{idx}", "-i", "{se}", "-r", "@RG\\tID:" + "x" * 256], "ID longer than 255"),
	("output-dir", ["-I", "{idx}", "-i", "{se}", "-o", "/nonexistent-dir/o.sam"], "Couldn't open output file"),
	("no-input", ["-I", "{idx}", "-i", "/nonexistent-dir/r.fq"], "Could not open file"),
]


@pytest.mark.parametrize("name,args,message", ERRORS, ids=[e[0] for e in ERRORS])
def test_argument_errors(gem3, se, pe, name, args, message):
	"""Invalid command lines exit 1 with the message, and write no SAM."""
	subst = {"{se}": se.paths[0], "{idx}": se.index, "{p1}": pe.paths[0], "{p2}": pe.paths[1]}
	argv = [subst.get(a, a) for a in args] + ["-t", "1"]
	res = mapper(gem3, *argv, timeout=60)
	assert res.returncode == 1 and not res.timed_out, res.describe()
	assert message in res.stderr_text, res.stderr_text[-2000:]
	assert res.stdout == b""


def test_rg_id_of_255_characters_is_accepted(se):
	"""A read-group ID of exactly 255 characters is accepted."""
	res = se.run(["-r", "@RG\\tID:" + "y" * 255], threads=1)
	assert res.sam().records[0].tags["RG"] == "y" * 255


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

# Every default printed by --help (except the known-bug ones) given
# explicitly, plus percentages written as such.
DEFAULTS = ["--clipping=uncalled", "-q", "offset-33", "--mapping-mode", "fast", "-E", "0.04",
	"-s", "1", "-e", "0.12", "--alignment-max-bandwidth", "0.20",
	"--alignment-global-min-identity", "80%", "--alignment-local", "if-unmapped",
	"--alignment-local-min-identity", "40", "--alignment-local-min-score", "20",
	"--alignment-force-full-swg=false", "--alignment-scaffolding-min-matching_length", "10",
	"--alignment-curation=true", "-A", "1", "-O", "6", "-X", "1", "--alignment-model", "gap-affine",
	"--mapq-model", "gem", "-M", "5", "-F", "SAM", "--sam-compact=true", "--sam-gem-compatible=false",
	"--map-format", "2", "--output-model", "4M,5c", "--input-model", "32M,c,4M",
	"--progress-step", "100000", "--bisulfite-read", "inferred",
	"--underconversion_sequence", "NC_001416.1", "--overconversion_sequence", "NC_001604.1",
	"--control_sequence", "NC_001422.1"]
PERCENT = ["-E", "4%", "-e", "12%", "--alignment-max-bandwidth", "20%"]
PE_DEFAULTS = ["-L", "10000", "--discordant-pair-search", "if-no-concordant"]


def test_documented_defaults_are_the_defaults(se, pe):
	"""Giving every documented default explicitly changes no output byte
	(mapper_parameters_set_defaults against the --help text)."""
	for ds, extra in ((se, []), (pe, PE_DEFAULTS)):
		ref = sam_digest(ds.run().stdout)
		assert sam_digest(ds.run(DEFAULTS + extra).stdout) == ref
		assert sam_digest(ds.run(PERCENT).stdout) == ref


def test_system_options_do_not_change_output(se, tmp_path):
	"""--max-memory and --tmp-folder are accepted and change no output."""
	ref = sam_digest(se.run().stdout)
	assert sam_digest(se.run(["--max-memory", "2GB"]).stdout) == ref
	assert sam_digest(se.run(["--tmp-folder", str(tmp_path)]).stdout) == ref


# ---------------------------------------------------------------------------
# Log output
# ---------------------------------------------------------------------------

def log_messages(res):
	out = []
	for line in res.stderr_text.split("\n"):
		if not line:
			continue
		m = LOG_LINE.match(line)
		assert m, line
		out.append(m.group(1))
	return out


def test_verbose_user_log(gem3, se, pe):
	"""--verbose=user logs the steps with a timestamp, and the total number
	of reads (single end) or pairs (paired end) processed."""
	for ds, label, n in ((se, "SE", len(se.reads or [])), (pe, "PE", len(pe.pairs or []))):
		res = ds.run(["--verbose=user"], threads=2)
		msgs = log_messages(res)
		opened = ["[Opening input file '%s']" % p for p in ds.paths]
		assert msgs[:len(opened)] == opened
		rest = msgs[len(opened):]
		assert rest[:3] == ["[Outputting to stdout]", "[Loading GEM index '%s']" % ds.index,
			"[%s::Mapping Sequences]" % label]
		assert rest[3] == "Total %d sequences processed... done" % n
		assert re.match(r"\A\[GEMMapper terminated successfully in \d+s\. \(\+\d+s\. loading\)\]\Z", rest[4])
		assert len(rest) == 5
	# -v without a value means user; '-v' then 'quiet' leaves 'quiet' unused.
	assert log_messages(se.run(["-v"], threads=2))[-1].startswith("[GEMMapper terminated")


def test_verbose_quiet_and_dev(se):
	"""quiet logs nothing; dev adds the index description."""
	assert se.run([], threads=1).stderr == b""
	dev = se.run(["--verbose=dev"], threads=1).stderr_text
	assert "[GEM]>Archive" in dev and "[GEM]>Locator" in dev
	ref = strip_variable_header(se.run([], threads=1).stdout)
	assert strip_variable_header(se.run(["--verbose=dev"], threads=1).stdout) == ref


@pytest.mark.parametrize("step", [25, 64])
def test_progress_step(se, step):
	"""--progress-step N logs '# k sequences processed' at each multiple of
	N (one mapping thread), then the total."""
	n = len(se.reads)
	msgs = log_messages(se.run(["--verbose=user", "--progress-step", str(step)], threads=1))
	ticks = [m for m in msgs if m.startswith("# ")]
	assert ticks == ["# %d sequences processed" % k for k in range(step, n, step)]
	assert "Total %d sequences processed... done" % n in msgs


def test_output_to_file_is_logged(gem3, se, tmp_path):
	res = gem3.run_map(se.index, inputs=se.paths, args=["--verbose=user", "-o", "x.sam"],
		cwd=str(tmp_path), threads=1)
	assert "[Outputting to 'x.sam']" in log_messages(res)


@pytest.mark.parametrize("reduce", ["sum", "min", "max", "mean", "sample"])
def test_profile(se, tmp_path, reduce):
	"""--profile=<reduce> prints the system information (date, working
	directory, host, user, CPUs) on stderr and changes no output."""
	res = se.run(["--profile=%s" % reduce], threads=2)
	text = res.stderr_text
	assert "[GEM]>System.Info\n" in text
	assert "  => CPUs %d\n" % os.sysconf("SC_NPROCESSORS_ONLN") in text
	assert "  => Username %s\n" % pwd.getpwuid(os.geteuid()).pw_name in text
	assert re.search(r"  => Date \d{4}/\d{1,2}/\d{1,2} \d\d:\d\d:\d\d CET\n", text)
	assert sam_digest(res.stdout) == sam_digest(se.run(threads=2).stdout)

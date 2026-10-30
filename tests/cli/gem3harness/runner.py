# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Running gem-indexer and gem-mapper under test.

Every process runs with ``LC_ALL=C``, core dumps off (``RLIMIT_CORE=1``,
which also keeps the kernel from piping a core to apport), a timeout, and a
private temporary working directory that is removed afterwards.
"""

import hashlib
import os
import resource
import shutil
import signal
import subprocess
import tempfile
import time

from .sam import parse_sam

DEFAULT_TIMEOUT = 120


def default_threads():
	"""-t for tests that do not pin it: GEM3_TEST_THREADS, default 1."""
	return int(os.environ.get("GEM3_TEST_THREADS", "1"))


def _limit_core():
	# RLIMIT_CORE=1 rather than 0: with a piped core_pattern (apport) the
	# kernel ignores a limit of 0 but skips the dump at 1.
	try:
		soft, hard = resource.getrlimit(resource.RLIMIT_CORE)
		if hard == resource.RLIM_INFINITY or hard >= 1:
			resource.setrlimit(resource.RLIMIT_CORE, (1, hard))
		else:
			resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
	except (ValueError, OSError):
		pass


def _env(extra=None):
	env = dict(os.environ)
	env["LC_ALL"] = "C"
	env.pop("LANG", None)
	if extra:
		env.update(extra)
	return env


class RunResult(object):
	"""``returncode``, ``stdout`` and ``stderr`` (bytes), ``seconds``,
	``files`` (name -> bytes of every file the process left in its working
	directory, when collected) and ``argv``. A negative returncode is the
	signal that killed the process; ``timed_out`` is set on a timeout."""

	def __init__(self, argv, returncode, stdout, stderr, seconds, files, timed_out=False):
		self.argv = argv
		self.returncode = returncode
		self.stdout = stdout
		self.stderr = stderr
		self.seconds = seconds
		self.files = files
		self.timed_out = timed_out

	@property
	def stderr_text(self):
		return self.stderr.decode("utf-8", "replace")

	@property
	def stdout_text(self):
		return self.stdout.decode("utf-8", "replace")

	@property
	def crashed(self):
		return self.returncode < 0 or self.timed_out

	def sam(self, name=None):
		"""Parse stdout (or the named output file) as SAM."""
		data = self.stdout if name is None else self.files[name]
		return parse_sam(data)

	def describe(self):
		return "argv: %s\nexit: %s%s\nstderr:\n%s" % (" ".join(self.argv), self.returncode,
			" (timed out)" if self.timed_out else "", self.stderr_text[-3000:])


def run_process(argv, cwd=None, stdin=None, timeout=DEFAULT_TIMEOUT, env=None, collect=True):
	"""Run ``argv`` in ``cwd`` (default: a fresh temporary directory that is
	deleted afterwards). ``stdin`` is bytes or None. Returns :class:`RunResult`."""
	own_dir = cwd is None
	if own_dir:
		cwd = tempfile.mkdtemp(prefix="gem3-run-")
	t0 = time.time()
	timed_out = False
	try:
		proc = subprocess.Popen(argv, cwd=cwd, stdin=subprocess.PIPE if stdin is not None else subprocess.DEVNULL,
			stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=_env(env), preexec_fn=_limit_core,
			start_new_session=True)
		try:
			out, err = proc.communicate(stdin, timeout=timeout)
		except subprocess.TimeoutExpired:
			timed_out = True
			try:
				os.killpg(proc.pid, signal.SIGKILL)
			except OSError:
				proc.kill()
			out, err = proc.communicate()
		files = {}
		if collect and own_dir:
			for name in sorted(os.listdir(cwd)):
				p = os.path.join(cwd, name)
				if os.path.isfile(p):
					with open(p, "rb") as fh:
						files[name] = fh.read()
		return RunResult(list(argv), proc.returncode, out, err, time.time() - t0, files, timed_out)
	finally:
		if own_dir:
			shutil.rmtree(cwd, ignore_errors=True)


class Gem3(object):
	"""Runs one build's tools. ``index()`` caches indexes per genome and
	options for the life of the object (one pytest session)."""

	def __init__(self, build, work_dir):
		self.build = build
		self.work_dir = work_dir
		self.index_seconds = {}
		self._indexes = {}
		os.makedirs(work_dir, exist_ok=True)

	@property
	def mapper(self):
		return self.build.tool("gem-mapper")

	@property
	def indexer(self):
		return self.build.tool("gem-indexer")

	def build_index(self, fasta, bisulfite=False, threads=2, timeout=300, extra_args=()):
		"""Index ``fasta`` (a path) with the gem-indexer under test; returns the
		.gem path. Cached by the FASTA's bytes and the options.

		gem-indexer writes ``<name>.sa.tmp`` into its working directory, so it
		runs with cwd set to the index's own directory.
		"""
		with open(fasta, "rb") as fh:
			digest = hashlib.sha256(fh.read()).hexdigest()[:16]
		key = (digest, bool(bisulfite), tuple(extra_args))
		if key in self._indexes:
			return self._indexes[key]
		name = "idx-%s%s" % (digest, "-bs" if bisulfite else "")
		if extra_args:
			name += "-" + hashlib.sha256(repr(extra_args).encode()).hexdigest()[:8]
		out_dir = os.path.join(self.work_dir, name)
		os.makedirs(out_dir, exist_ok=True)
		prefix = os.path.join(out_dir, "index")
		argv = [self.indexer, "-i", os.path.abspath(fasta), "-o", prefix, "-t", str(threads)]
		if bisulfite:
			argv.append("-b")
		argv += list(extra_args)
		res = run_process(argv, cwd=out_dir, timeout=timeout, collect=False)
		if res.returncode != 0 or not os.path.exists(prefix + ".gem"):
			raise RuntimeError("gem-indexer failed\n" + res.describe())
		self.index_seconds[key] = res.seconds
		self._indexes[key] = prefix + ".gem"
		return prefix + ".gem"

	def index(self, genome, bisulfite=False, **kw):
		"""Write a :class:`Genome` to FASTA in the work dir and index it."""
		fa = os.path.join(self.work_dir, "genome-%s.fa" % genome.digest()[:16])
		if not os.path.exists(fa):
			genome.write_fasta(fa + ".tmp")
			os.rename(fa + ".tmp", fa)
		return self.build_index(fa, bisulfite=bisulfite, **kw)

	def run_map(self, index, args=(), inputs=(), stdin=None, threads=None, timeout=DEFAULT_TIMEOUT,
			check=True, quiet=True, files=None, cwd=None):
		"""Run gem-mapper on ``index``; returns :class:`RunResult`.

		``args`` are extra arguments. ``inputs`` are paths passed with -i (one)
		or -1/-2 (two). ``files`` maps names to bytes written into the private
		working directory first, so arguments can name them relatively (e.g.
		``files={"r.fq": data}, args=["-i", "r.fq"]``). ``threads`` defaults to
		GEM3_TEST_THREADS; pass ``threads=False`` to leave -t out. With
		``check`` a non-zero exit raises AssertionError with stderr.
		"""
		argv = [self.mapper, "-I", os.path.abspath(index)]
		inputs = list(inputs)
		if len(inputs) == 1:
			argv += ["-i", os.path.abspath(inputs[0])]
		elif len(inputs) == 2:
			argv += ["-1", os.path.abspath(inputs[0]), "-2", os.path.abspath(inputs[1])]
		elif inputs:
			raise ValueError("inputs takes one or two paths")
		if threads is None:
			threads = default_threads()
		if threads is not False:
			argv += ["-t", str(threads)]
		# -v takes an optional argument, which getopt only binds when attached
		# ("-v quiet" sets the default verbosity and leaves "quiet" as an
		# ignored positional argument).
		if quiet and not any(str(a).startswith("-v") or str(a).startswith("--verbose") for a in args):
			argv += ["--verbose=quiet"]
		argv += [str(a) for a in args]
		own = cwd is None
		if files or not own:
			d = cwd or tempfile.mkdtemp(prefix="gem3-map-")
			try:
				for name, data in (files or {}).items():
					with open(os.path.join(d, name), "wb") as fh:
						fh.write(data)
				res = run_process(argv, cwd=d, stdin=stdin, timeout=timeout, collect=False)
				res.files = {}
				for name in sorted(os.listdir(d)):
					p = os.path.join(d, name)
					if os.path.isfile(p) and name not in (files or {}):
						with open(p, "rb") as fh:
							res.files[name] = fh.read()
			finally:
				if own:
					shutil.rmtree(d, ignore_errors=True)
		else:
			res = run_process(argv, stdin=stdin, timeout=timeout)
		if check and (res.returncode != 0 or res.timed_out):
			raise AssertionError("gem-mapper failed\n" + res.describe())
		return res

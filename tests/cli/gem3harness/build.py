# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""Out-of-tree builds of gem3-mapper and of the C unit-test drivers.

``build_gem3()`` copies the source tree (everything but ``.git``, ``tests``
and build products) into ``<cache>/build-<key>``, where the key hashes the
sources, the compiler and the make arguments, then runs ``./configure
--disable-cuda`` and ``make`` there. Nothing is written to the source tree.
A file lock serializes concurrent builders; a finished build is reused.

The tree's Makefiles race under ``make -jN`` (the tools link while the
objects are still compiling), so the build compiles with ``make -k -jN``,
deletes whatever was linked, and links serially with a second ``make``.
This produces the same binaries as a serial ``make``.

Environment:
  GEM3_BIN_DIR       use the binaries in this directory and build nothing
                     (C drivers are skipped unless GEM3_BUILD_DIR is set too)
  GEM3_BUILD_DIR     an existing configured-and-built gem3 tree whose
                     ``build/*.o`` the C drivers link against
  GEM3_TEST_CACHE    cache directory (default $XDG_CACHE_HOME/gem3-mapper-tests
                     or ~/.cache/gem3-mapper-tests)
  GEM3_CC            compiler passed to configure as CC (default: configure's choice)
  GEM3_MAKE_TARGET   release (default), debug, devel or profile
  GEM3_MAKE_ARGS     extra make arguments, e.g. "FLAGS_OPT=-O3"
  GEM3_BUILD_JOBS    parallel compile jobs (default 4)
"""

import fcntl
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

TESTS_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REPO_ROOT = os.path.dirname(TESTS_DIR)

# Top-level names never copied or hashed: VCS data, the test suite itself and
# what configure/make write.
_EXCLUDE_TOP = {".git", "tests", "bin", "build", "lib", "Makefile.mk", "config.status",
	"config.log", "configure.log", ".project", ".cproject", ".settings"}
_EXCLUDE_ANY = {"__pycache__", ".pytest_cache"}
TOOLS = ("gem-mapper", "gem-indexer", "gem-retriever")


class BuildError(RuntimeError):
	pass


def cache_root():
	root = os.environ.get("GEM3_TEST_CACHE")
	if not root:
		base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
		root = os.path.join(base, "gem3-mapper-tests")
	os.makedirs(root, exist_ok=True)
	return root


def _source_entries(root):
	"""Sorted (relpath, abspath) of every source file to copy."""
	out = []
	for dirpath, dirnames, filenames in os.walk(root):
		rel_dir = os.path.relpath(dirpath, root)
		if rel_dir == ".":
			dirnames[:] = [d for d in dirnames if d not in _EXCLUDE_TOP and d not in _EXCLUDE_ANY]
			filenames = [f for f in filenames if f not in _EXCLUDE_TOP]
		else:
			dirnames[:] = [d for d in dirnames if d not in _EXCLUDE_ANY]
		# include/resources is the symlink `make setup` creates in a tree
		# that was built in place.
		if rel_dir == "include":
			dirnames[:] = [d for d in dirnames if d != "resources"]
			filenames = [f for f in filenames if f != "resources"]
		# resources/build is created by the resources Makefile.
		if rel_dir == "resources":
			dirnames[:] = [d for d in dirnames if d != "build"]
		for d in list(dirnames):
			if os.path.islink(os.path.join(dirpath, d)):
				dirnames.remove(d)
				filenames.append(d)
		for f in filenames:
			p = os.path.join(dirpath, f)
			out.append((os.path.normpath(os.path.join(rel_dir, f)), p))
	out.sort()
	return out


def source_digest(root=REPO_ROOT):
	"""sha256 over the relative path, mode and content of every source file."""
	h = hashlib.sha256()
	for rel, path in _source_entries(root):
		h.update(rel.encode("utf-8") + b"\0")
		if os.path.islink(path):
			h.update(b"L" + os.readlink(path).encode("utf-8") + b"\0")
			continue
		st = os.stat(path)
		h.update(b"X" if st.st_mode & 0o111 else b"F")
		with open(path, "rb") as fh:
			h.update(hashlib.sha256(fh.read()).digest())
	return h.hexdigest()


def tree_snapshot(root=REPO_ROOT):
	"""(relpath -> (size, mtime_ns)) of the source tree, tests/ excluded, to
	show a build left it untouched."""
	snap = {}
	for dirpath, dirnames, filenames in os.walk(root):
		rel_dir = os.path.relpath(dirpath, root)
		if rel_dir == ".":
			dirnames[:] = [d for d in dirnames if d not in (".git", "tests")]
		for f in filenames + [d for d in dirnames if os.path.islink(os.path.join(dirpath, d))]:
			p = os.path.join(dirpath, f)
			st = os.lstat(p)
			snap[os.path.normpath(os.path.join(rel_dir, f))] = (st.st_size, st.st_mtime_ns)
	return snap


def _copy_sources(root, dest):
	for rel, path in _source_entries(root):
		target = os.path.join(dest, rel)
		os.makedirs(os.path.dirname(target), exist_ok=True)
		if os.path.islink(path):
			os.symlink(os.readlink(path), target)
		else:
			shutil.copy2(path, target)


def _compiler_id(cc):
	try:
		out = subprocess.run([cc, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
			universal_newlines=True, timeout=60).stdout
	except (OSError, subprocess.SubprocessError) as exc:
		raise BuildError("cannot run compiler %r: %s" % (cc, exc))
	return out.strip()


class Gem3Build(object):
	"""A finished build: ``bin_dir`` with the tools, ``root`` with ``build/*.o``.

	``cc``, ``compile_flags`` and ``link_flags`` are what the tree's Makefile
	used to link gem-mapper (parsed from the make log), for the C drivers.
	"""

	def __init__(self, root, bin_dir, info):
		self.root = root
		self.bin_dir = bin_dir
		self.info = info

	@property
	def cc(self):
		return self.info.get("cc")

	@property
	def compile_flags(self):
		return self.info.get("compile_flags", [])

	@property
	def link_flags(self):
		return self.info.get("link_flags", [])

	@property
	def has_objects(self):
		return self.root is not None and os.path.isdir(os.path.join(self.root, "build"))

	def tool(self, name):
		return os.path.join(self.bin_dir, name)

	def describe(self):
		if self.root is None:
			return "binaries from GEM3_BIN_DIR=%s" % self.bin_dir
		return "%s %s (%s)" % (self.info.get("target"), " ".join(self.info.get("make_args", [])),
			self.root)


def _parse_link_line(log_text, root):
	"""Find the command that linked gem-mapper and split it into compiler,
	compile flags and link flags (paths made absolute)."""
	for line in log_text.splitlines():
		toks = shlex.split(line) if "gem-mapper.c" in line and " -o " in line else None
		if not toks or "gem-mapper.c" not in toks:
			continue
		cc = toks[0]
		src_idx = toks.index("gem-mapper.c")
		compile_flags = [t for t in toks[1:src_idx] if not t.startswith("-DGEM_VERSION")]
		o_idx = toks.index("-o")
		tail = toks[o_idx + 2:]
		link = []
		for t in tail:
			if t.startswith("-I") or t.startswith("-L"):
				continue
			link.append(t)
		return cc, compile_flags, link
	raise BuildError("could not find the gem-mapper link command in the make log")


def build_gem3(cc=None, target=None, make_args=None, jobs=None, root=REPO_ROOT, log=None):
	"""Build (or reuse) gem3-mapper out of tree; returns a :class:`Gem3Build`."""
	cc = cc if cc is not None else os.environ.get("GEM3_CC") or None
	target = target or os.environ.get("GEM3_MAKE_TARGET") or "release"
	if make_args is None:
		make_args = shlex.split(os.environ.get("GEM3_MAKE_ARGS", ""))
	jobs = int(jobs or os.environ.get("GEM3_BUILD_JOBS") or 4)
	cc_for_id = cc or "gcc"
	key_src = json.dumps({
		"sources": source_digest(root),
		"cc": cc,
		"cc_version": _compiler_id(cc_for_id),
		"target": target,
		"make_args": list(make_args),
		"layout": 1,
	}, sort_keys=True)
	key = hashlib.sha256(key_src.encode("utf-8")).hexdigest()[:20]
	croot = cache_root()
	bdir = os.path.join(croot, "build-" + key)
	lock_path = bdir + ".lock"
	with open(lock_path, "w") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		info_path = os.path.join(bdir, ".gem3-test-build.json")
		if os.path.exists(info_path):
			with open(info_path) as fh:
				info = json.load(fh)
			info["reused"] = True
			return Gem3Build(bdir, os.path.join(bdir, "bin"), info)
		if os.path.exists(bdir):
			shutil.rmtree(bdir)
		os.makedirs(bdir)
		t0 = time.time()
		_copy_sources(root, bdir)
		env = dict(os.environ)
		env["LC_ALL"] = "C"
		if cc:
			env["CC"] = cc
		logs = []

		def run(cmd, check=True):
			proc = subprocess.run(cmd, cwd=bdir, env=env, stdout=subprocess.PIPE,
				stderr=subprocess.STDOUT, universal_newlines=True)
			logs.append("$ " + " ".join(cmd) + "\n" + proc.stdout)
			if check and proc.returncode != 0:
				with open(os.path.join(bdir, "build.log"), "w") as fh:
					fh.write("\n".join(logs))
				raise BuildError("%s failed (exit %d); log in %s/build.log\n%s"
					% (" ".join(cmd), proc.returncode, bdir, proc.stdout[-4000:]))
			return proc

		run(["./configure", "--disable-cuda"])
		run(["make", "-k", "-j%d" % jobs, target] + list(make_args), check=False)
		for d in ("bin", "lib"):
			p = os.path.join(bdir, d)
			if os.path.isdir(p):
				for f in os.listdir(p):
					os.remove(os.path.join(p, f))
		link = run(["make", target] + list(make_args))
		for tool in TOOLS:
			if not os.access(os.path.join(bdir, "bin", tool), os.X_OK):
				raise BuildError("make finished but bin/%s is missing" % tool)
		cc_used, compile_flags, link_flags = _parse_link_line(link.stdout, bdir)
		info = {
			"key": key,
			"cc": cc_used,
			"cc_requested": cc,
			"cc_version": _compiler_id(cc_used),
			"target": target,
			"make_args": list(make_args),
			"compile_flags": compile_flags,
			"link_flags": link_flags,
			"build_seconds": round(time.time() - t0, 2),
			"jobs": jobs,
		}
		with open(os.path.join(bdir, "build.log"), "w") as fh:
			fh.write("\n".join(logs))
		with open(info_path + ".tmp", "w") as fh:
			json.dump(info, fh, indent=1)
		os.rename(info_path + ".tmp", info_path)
		info["reused"] = False
		if log:
			log("built gem3-mapper in %.1f s at %s" % (info["build_seconds"], bdir))
		return Gem3Build(bdir, os.path.join(bdir, "bin"), info)


def gem3_from_env(log=None):
	"""The build under test: GEM3_BIN_DIR if set, else :func:`build_gem3`."""
	bin_dir = os.environ.get("GEM3_BIN_DIR")
	if bin_dir:
		bin_dir = os.path.abspath(bin_dir)
		for tool in ("gem-mapper", "gem-indexer"):
			if not os.access(os.path.join(bin_dir, tool), os.X_OK):
				raise BuildError("GEM3_BIN_DIR=%s has no executable %s" % (bin_dir, tool))
		root = os.environ.get("GEM3_BUILD_DIR")
		info = {"reused": True, "target": "external"}
		if root:
			root = os.path.abspath(root)
			info.update(_flags_from_tree(root))
		return Gem3Build(root, bin_dir, info)
	return build_gem3(log=log)


def _flags_from_tree(root):
	"""Ask an already-built tree's make how it would link gem-mapper."""
	# A dry run of the whole tree (-n -B) prints the gem-mapper link command
	# with the flags the Makefiles pass down; nothing is written.
	target = os.environ.get("GEM3_MAKE_TARGET") or "release"
	proc = subprocess.run(["make", "-n", "-B", "-C", root, target], stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT, universal_newlines=True)
	cc, compile_flags, link_flags = _parse_link_line(proc.stdout, root)
	return {"cc": cc, "compile_flags": compile_flags, "link_flags": link_flags}


# --------------------------------------------------------------------------
# C unit-test drivers
# --------------------------------------------------------------------------

C_DIR = os.path.join(TESTS_DIR, "c")


def c_driver_sources():
	"""tests/c/test_*.c, sorted."""
	return sorted(os.path.join(C_DIR, f) for f in os.listdir(C_DIR)
		if f.startswith("test_") and f.endswith(".c"))


def build_c_driver(gem3, source, log=None):
	"""Compile and link one driver against ``gem3``'s objects; returns its path.

	The driver is compiled with the flags the tree used for gem-mapper and
	linked with every object in ``build/`` (the tools' own ``main``s are not
	in there). Cached by the driver's source, the harness header and the build.
	"""
	if not gem3.has_objects or not gem3.cc:
		raise BuildError("no gem3 object directory to link C drivers against "
			"(set GEM3_BUILD_DIR when GEM3_BIN_DIR is set)")
	objs = sorted(os.path.join(gem3.root, "build", f)
		for f in os.listdir(os.path.join(gem3.root, "build")) if f.endswith(".o"))
	h = hashlib.sha256()
	for p in [source, os.path.join(C_DIR, "gem3_test.h")]:
		with open(p, "rb") as fh:
			h.update(hashlib.sha256(fh.read()).digest())
	h.update(json.dumps([gem3.root, gem3.info.get("key"), gem3.cc, gem3.compile_flags,
		gem3.link_flags]).encode("utf-8"))
	for o in objs:
		st = os.stat(o)
		h.update(("%s %d %d" % (o, st.st_size, st.st_mtime_ns)).encode("utf-8"))
	key = h.hexdigest()[:20]
	name = os.path.splitext(os.path.basename(source))[0]
	out_dir = os.path.join(cache_root(), "c-drivers")
	os.makedirs(out_dir, exist_ok=True)
	exe = os.path.join(out_dir, "%s-%s" % (name, key))
	with open(exe + ".lock", "w") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		if os.access(exe, os.X_OK):
			return exe
		cmd = ([gem3.cc] + list(gem3.compile_flags) + ["-I" + os.path.join(gem3.root, "include"),
			"-I" + C_DIR, source] + objs + ["-o", exe + ".tmp"] + list(gem3.link_flags))
		t0 = time.time()
		proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
			universal_newlines=True, env=dict(os.environ, LC_ALL="C"))
		if proc.returncode != 0:
			shown = [c for c in cmd if c not in objs]
			shown.insert(shown.index("-o"), "<%d objects from %s/build>" % (len(objs), gem3.root))
			raise BuildError("compiling %s failed:\n%s\n%s" % (name, " ".join(shown), proc.stdout[-6000:]))
		os.rename(exe + ".tmp", exe)
		if log:
			log("built %s in %.1f s" % (name, time.time() - t0))
	return exe


if __name__ == "__main__":
	# python3 tests/cli/gem3harness/build.py  -> build and print the build directory
	b = gem3_from_env(log=lambda m: print(m, file=sys.stderr))
	print(b.bin_dir)

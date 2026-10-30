# GEM-Mapper v3 (GEM3) test suite.
# Copyright (c) 2026 Jacob Schreiber <jacob.schreiber@umassmed.edu>
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of GEM-Mapper v3 (GEM3), distributed under the GNU General
# Public License version 3 or (at your option) any later version; see LICENSE.
"""A small SAM parser and independent recomputation of CIGAR, NM and MD.

Written from the SAM specification (SAMv1, sections 1.4 and 1.5), not from
gem3-mapper's writer. ``parse_sam`` returns a :class:`Sam` with the header
and :class:`SamRecord` objects whose optional tags are typed (``i`` -> int,
``f`` -> float, ``A``/``Z``/``H`` -> str, ``B`` -> (subtype, list)).
"""

import re

# SAM FLAG bits
PAIRED = 0x1
PROPER_PAIR = 0x2
UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
REVERSE = 0x10
MATE_REVERSE = 0x20
READ1 = 0x40
READ2 = 0x80
SECONDARY = 0x100
QCFAIL = 0x200
DUPLICATE = 0x400
SUPPLEMENTARY = 0x800

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_MD_RE = re.compile(r"\A[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*\Z")
_MD_TOKEN_RE = re.compile(r"(\d+)|(\^[A-Z]+)|([A-Z])")

CONSUMES_QUERY = set("MIS=X")
CONSUMES_REF = set("MDN=X")


class SamFormatError(ValueError):
	pass


def parse_cigar(cigar):
	"""'10M2I5M' -> [('M', 10), ('I', 2), ('M', 5)]; '*' -> []."""
	if cigar == "*":
		return []
	ops = [(op, int(n)) for n, op in _CIGAR_RE.findall(cigar)]
	if "".join("%d%s" % (n, op) for op, n in ops) != cigar:
		raise SamFormatError("malformed CIGAR %r" % cigar)
	return ops


def cigar_to_string(ops):
	return "".join("%d%s" % (n, op) for op, n in ops) if ops else "*"


def cigar_query_length(ops):
	"""Read bases the CIGAR accounts for (M, I, S, =, X)."""
	return sum(n for op, n in ops if op in CONSUMES_QUERY)


def cigar_ref_span(ops):
	"""Reference bases the CIGAR covers (M, D, N, =, X)."""
	return sum(n for op, n in ops if op in CONSUMES_REF)


def cigar_is_valid(ops, seq_length=None):
	"""SAM rules: H only at the ends, S only next to the ends or an H, no
	zero lengths, at least one reference-consuming op, and the query length
	equal to SEQ's length when given."""
	if not ops:
		return False
	if any(n <= 0 for _, n in ops):
		return False
	inner = list(ops)
	if inner and inner[0][0] == "H":
		inner = inner[1:]
	if inner and inner[-1][0] == "H":
		inner = inner[:-1]
	if any(op == "H" for op, _ in inner):
		return False
	core = list(inner)
	if core and core[0][0] == "S":
		core = core[1:]
	if core and core[-1][0] == "S":
		core = core[:-1]
	if any(op == "S" for op, _ in core):
		return False
	if not any(op in CONSUMES_REF for op, _ in core):
		return False
	if seq_length is not None and cigar_query_length(ops) != seq_length:
		return False
	return True


def md_nm_from_alignment(ref, pos, ops, seq):
	"""Recompute (MD, NM) from the reference, POS (0-based), CIGAR and SEQ.

	Comparison is case-insensitive (the MD string records the upper-case
	reference base). A read N against a reference base counts as a mismatch,
	and so does any read base against a reference N. Soft clips are skipped;
	``N`` (skipped region) and ``P`` add nothing; ``=``/``X`` are compared
	like ``M``. Returns (md_string, nm).
	"""
	md = []
	run = 0
	nm = 0
	r = pos
	q = 0
	for op, n in ops:
		if op in ("M", "=", "X"):
			for k in range(n):
				rb = ref[r + k].upper()
				qb = seq[q + k].upper()
				if rb == qb and rb != "N":
					run += 1
				else:
					md.append(str(run))
					md.append(rb)
					run = 0
					nm += 1
			r += n
			q += n
		elif op == "I":
			nm += n
			q += n
		elif op == "D":
			md.append(str(run))
			md.append("^" + ref[r:r + n].upper())
			run = 0
			nm += n
			r += n
		elif op == "N":
			r += n
		elif op == "S":
			q += n
		elif op in ("H", "P"):
			pass
		else:
			raise SamFormatError("unknown CIGAR op %r" % op)
	md.append(str(run))
	return "".join(md), nm


def bisulfite_nm(ref, pos, ops, seq, conversion):
	"""Edit distance that ignores bisulfite conversions.

	``conversion`` is 'C' (C-to-T on the forward reference strand: a read T
	over a reference C is not counted) or 'G' (G-to-A: a read A over a
	reference G is not counted). Every other mismatch, inserted and deleted
	base counts once.
	"""
	allowed = {"C": ("C", "T"), "G": ("G", "A")}[conversion]
	nm = 0
	r = pos
	q = 0
	for op, n in ops:
		if op in ("M", "=", "X"):
			for k in range(n):
				rb = ref[r + k].upper()
				qb = seq[q + k].upper()
				if rb == qb and rb != "N":
					continue
				if (rb, qb) == allowed:
					continue
				nm += 1
			r += n
			q += n
		elif op == "I":
			nm += n
			q += n
		elif op == "D":
			nm += n
			r += n
		elif op == "N":
			r += n
		elif op == "S":
			q += n
	return nm


def md_is_wellformed(md):
	return bool(_MD_RE.match(md))


def md_ref_span(md):
	"""Reference bases an MD string covers (matches + mismatches + deletions)."""
	span = 0
	for num, dele, mis in _MD_TOKEN_RE.findall(md):
		if num:
			span += int(num)
		elif dele:
			span += len(dele) - 1
		else:
			span += 1
	return span


def md_equal(a, b):
	"""MD strings equal up to the optional zero-length match runs SAMv1 allows."""
	def norm(s):
		toks = []
		for num, dele, mis in _MD_TOKEN_RE.findall(s):
			if num:
				toks.append(int(num))
			else:
				toks.append(dele or mis)
		return toks
	return norm(a) == norm(b)


class SamRecord(object):
	"""One alignment line. ``pos`` and ``pnext`` are 1-based as in the file."""

	__slots__ = ("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext",
		"tlen", "seq", "qual", "tags", "tag_types", "line")

	def __init__(self, line):
		fields = line.rstrip("\n").split("\t")
		if len(fields) < 11:
			raise SamFormatError("SAM record with %d fields: %r" % (len(fields), line))
		self.line = line.rstrip("\n")
		self.qname = fields[0]
		self.flag = int(fields[1])
		self.rname = fields[2]
		self.pos = int(fields[3])
		self.mapq = int(fields[4])
		self.cigar = fields[5]
		self.rnext = fields[6]
		self.pnext = int(fields[7])
		self.tlen = int(fields[8])
		self.seq = fields[9]
		self.qual = fields[10]
		self.tags = {}
		self.tag_types = {}
		for field in fields[11:]:
			if len(field) < 5 or field[2] != ":" or field[4] != ":":
				raise SamFormatError("malformed optional field %r" % field)
			tag, typ, value = field[:2], field[3], field[5:]
			if tag in self.tags:
				raise SamFormatError("duplicate tag %s in %r" % (tag, line))
			self.tags[tag] = _typed(typ, value)
			self.tag_types[tag] = typ

	# Flag helpers
	def has(self, bit):
		return bool(self.flag & bit)

	@property
	def is_unmapped(self):
		return self.has(UNMAPPED)

	@property
	def is_reverse(self):
		return self.has(REVERSE)

	@property
	def is_secondary(self):
		return self.has(SECONDARY)

	@property
	def is_supplementary(self):
		return self.has(SUPPLEMENTARY)

	@property
	def is_primary(self):
		return not (self.flag & (SECONDARY | SUPPLEMENTARY))

	@property
	def strand(self):
		return "-" if self.is_reverse else "+"

	@property
	def cigar_ops(self):
		return parse_cigar(self.cigar)

	@property
	def ref_start(self):
		"""0-based leftmost reference position."""
		return self.pos - 1

	@property
	def ref_end(self):
		"""0-based exclusive end, from the CIGAR."""
		return self.pos - 1 + cigar_ref_span(self.cigar_ops)

	def tag(self, name, default=None):
		return self.tags.get(name, default)

	def __repr__(self):
		return "SamRecord(%s)" % self.line[:200]


def _typed(typ, value):
	if typ == "i":
		return int(value)
	if typ == "f":
		return float(value)
	if typ in ("A", "Z", "H"):
		return value
	if typ == "B":
		parts = value.split(",")
		sub = parts[0]
		conv = float if sub == "f" else int
		return (sub, [conv(x) for x in parts[1:]])
	raise SamFormatError("unknown tag type %r" % typ)


class Sam(object):
	"""Header lines and records of a SAM file."""

	def __init__(self, header_lines, records):
		self.header_lines = header_lines
		self.records = records

	@property
	def sq(self):
		"""[(name, length)] from @SQ lines, in order."""
		out = []
		for line in self.header_lines:
			if line.startswith("@SQ\t"):
				f = header_fields(line)
				out.append((f["SN"], int(f["LN"])))
		return out

	def header(self, record_type):
		"""All header lines of a type (e.g. '@PG') as dicts of their fields."""
		return [header_fields(l) for l in self.header_lines if l.startswith(record_type + "\t")]

	def by_name(self):
		"""qname -> list of records in file order."""
		out = {}
		for r in self.records:
			out.setdefault(r.qname, []).append(r)
		return out

	def primary(self):
		return [r for r in self.records if r.is_primary]


def header_fields(line):
	fields = line.rstrip("\n").split("\t")
	out = {"_type": fields[0]}
	for f in fields[1:]:
		if len(f) >= 3 and f[2] == ":":
			out[f[:2]] = f[3:]
	return out


def parse_sam(text):
	"""Parse SAM text (str or bytes)."""
	if isinstance(text, bytes):
		text = text.decode("ascii")
	header = []
	records = []
	for line in text.split("\n"):
		if not line:
			continue
		if line.startswith("@"):
			if records:
				raise SamFormatError("header line after records: %r" % line)
			header.append(line)
		else:
			records.append(SamRecord(line))
	return Sam(header, records)


# Header lines whose content depends on the run, not on the mapping: @PG
# (VN is the build's version string, CL the argv after getopt permuted it and
# the comma-separated option values were split in place) and @CO (TM time, WD
# working directory, HN host, UN user).
VARIABLE_HEADER_TYPES = ("@PG", "@CO")


def strip_variable_header(data):
	"""SAM bytes (or str) without the @PG and @CO lines, for byte comparisons
	and pinned digests."""
	text = data.decode("ascii") if isinstance(data, bytes) else data
	kept = [l for l in text.split("\n") if not l.startswith(VARIABLE_HEADER_TYPES)]
	out = "\n".join(kept)
	return out.encode("ascii") if isinstance(data, bytes) else out


def sam_digest(data):
	"""sha256 of SAM output with the variable header lines removed."""
	import hashlib
	stripped = strip_variable_header(data)
	if not isinstance(stripped, bytes):
		stripped = stripped.encode("ascii")
	return hashlib.sha256(stripped).hexdigest()


def check_record_against_reference(rec, ref_seq):
	"""Problems with a mapped record's CIGAR, NM and MD against ``ref_seq``.

	Returns a list of strings (empty when consistent): CIGAR valid and its
	query length equal to len(SEQ); the alignment inside the contig; NM and
	MD, when present, equal to the values recomputed from the reference.

	When SEQ is ``*`` (secondary records), NM and MD cannot be recomputed:
	the CIGAR must be valid, the alignment inside the contig, MD well formed
	with the CIGAR's reference span and its deleted bases equal to the
	reference, and NM equal to MD's mismatches plus the CIGAR's indel bases.
	"""
	problems = []
	ops = rec.cigar_ops
	no_seq = rec.seq == "*"
	if no_seq and not cigar_is_valid(ops):
		problems.append("CIGAR %s invalid" % rec.cigar)
		return problems
	if not no_seq and not cigar_is_valid(ops, len(rec.seq)):
		problems.append("CIGAR %s invalid for a %d bp SEQ" % (rec.cigar, len(rec.seq)))
		return problems
	if rec.ref_start < 0 or rec.ref_end > len(ref_seq):
		problems.append("alignment %d-%d outside contig of length %d"
			% (rec.ref_start, rec.ref_end, len(ref_seq)))
		return problems
	if no_seq:
		return problems + _md_nm_without_seq(rec, ref_seq, ops)
	md, nm = md_nm_from_alignment(ref_seq, rec.ref_start, ops, rec.seq)
	if "NM" in rec.tags and rec.tags["NM"] != nm:
		problems.append("NM:%d but recomputed %d" % (rec.tags["NM"], nm))
	if "MD" in rec.tags:
		if not md_is_wellformed(rec.tags["MD"]):
			problems.append("MD:%s malformed" % rec.tags["MD"])
		elif not md_equal(rec.tags["MD"], md):
			problems.append("MD:%s but recomputed %s" % (rec.tags["MD"], md))
	return problems


def _md_nm_without_seq(rec, ref_seq, ops):
	"""The NM/MD checks possible on a record whose SEQ is '*'."""
	problems = []
	indels = sum(n for op, n in ops if op in ("I", "D"))
	if "MD" not in rec.tags:
		if "NM" in rec.tags and rec.tags["NM"] < indels:
			problems.append("NM:%d below the CIGAR's %d indel bases" % (rec.tags["NM"], indels))
		return problems
	md = rec.tags["MD"]
	if not md_is_wellformed(md):
		return ["MD:%s malformed" % md]
	span = sum(n for op, n in ops if op in ("M", "=", "X", "D"))
	if md_ref_span(md) != span:
		problems.append("MD:%s spans %d reference bases, CIGAR %s %d"
			% (md, md_ref_span(md), rec.cigar, span))
		return problems
	# Walk MD along the CIGAR's reference positions: mismatches must sit on
	# M/=/X columns, deleted bases on D columns and equal the reference.
	cols = []
	r = rec.ref_start
	for op, n in ops:
		if op in ("M", "=", "X", "D"):
			cols.extend((op, r + k) for k in range(n))
			r += n
		elif op == "N":
			r += n
	i = 0
	mismatches = 0
	for num, dele, mis in _MD_TOKEN_RE.findall(md):
		if num:
			i += int(num)
		elif dele:
			for k, base in enumerate(dele[1:]):
				op, p = cols[i + k]
				if op != "D" or ref_seq[p].upper() != base:
					problems.append("MD:%s deletion base %s at %d is not a deleted reference base" % (md, base, p))
			i += len(dele) - 1
		else:
			op, p = cols[i]
			if op == "D" or ref_seq[p].upper() != mis:
				problems.append("MD:%s mismatch base %s at %d is not the reference base" % (md, mis, p))
			mismatches += 1
			i += 1
	if "NM" in rec.tags and rec.tags["NM"] != mismatches + indels:
		problems.append("NM:%d but MD's %d mismatches plus %d indel bases"
			% (rec.tags["NM"], mismatches, indels))
	return problems

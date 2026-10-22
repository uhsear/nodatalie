#!/usr/bin/env python
"""Gate on stored versus recomputed raster statistics.

A raster carries its minimum, maximum, mean and standard deviation as stored
numbers, in a PAM sidecar (name.tif.aux.xml) or in the file itself. Tools that
draw a legend, stretch a colour ramp or sanity-check an import read those stored
numbers. They do not scan the cells. Change the NoData value afterwards, or
overwrite the cells, and the stored numbers keep describing a raster that no
longer exists.

This tool reads the cells, recomputes the statistics with the NoData cells
masked, and compares them with the stored numbers. A mismatch is a finding. It
also names the case where the stored numbers counted the cells that the
declared NoData value now masks, and the case where the NoData value sits
inside the real data range and so masks real cells.

    python nodatalie.py --self-test
    python nodatalie.py rasters/
    python nodatalie.py dem.tif --json
    python nodatalie.py rasters/ --apply

Reads Esri ASCII grids (.asc) and single band GeoTIFFs (.tif, .tiff) with
strips or tiles, uncompressed, LZW or deflate, in either byte order. Standard
library only. Nothing is written without --apply, and --apply writes only the
PAM sidecar, never the raster.

Exit codes: 0 clean, 1 a finding, 2 a raster could not be read or a path is
missing, 64 usage error.
"""

from __future__ import print_function

import argparse
import array
import ast
import contextlib
import io
import itertools
import json
import math
import os
import shutil
import struct
import sys
import tempfile
import xml.etree.ElementTree as ET
import zlib

# =============================================================================
# CONFIGURATION. Deliberately not flags. Change here, not at the call site.
# =============================================================================

# Relative tolerance when a stored number is compared with a recomputed one.
# Stored numbers are text. GDAL writes 14 significant digits, so 1e-6 is loose
# enough for other writers and far tighter than any real change.
DEFAULT_TOL = 1e-6

# The stored valid-percent is rounded to four significant digits by GDAL, so
# it can only be compared to about this many percentage points.
VALID_PCT_TOL = 0.01

# Largest read from an uncompressed strip or an ASCII grid in one go.
BLOCK_BYTES = 1 << 20

# =============================================================================
# End of CONFIGURATION.
# =============================================================================

OK = "OK"
WARN = "WARN"
FAIL = "FAIL"
UNREADABLE = "UNREADABLE"

K_MIN = "STATISTICS_MINIMUM"
K_MAX = "STATISTICS_MAXIMUM"
K_MEAN = "STATISTICS_MEAN"
K_STD = "STATISTICS_STDDEV"
K_VALID = "STATISTICS_VALID_PERCENT"
K_APPROX = "STATISTICS_APPROXIMATE"
STAT_KEYS = ((K_MIN, "min"), (K_MAX, "max"), (K_MEAN, "mean"), (K_STD, "std"))

RASTER_EXT = (".tif", ".tiff", ".asc")
BAK_SUFFIX = ".nodatalie.bak"

ASC_KEYS = (b"ncols", b"nrows", b"xllcorner", b"yllcorner", b"xllcenter",
            b"yllcenter", b"cellsize", b"dx", b"dy", b"nodata_value")


class UnreadableRaster(Exception):
    """The file is not a raster this tool can read, or is damaged."""


# ----------------------------------------------------------------- pure core

EMPTY = (0, 0.0, 0.0, None, None)


def acc_of(values):
    """Accumulator (count, mean, sum of squared deviations, min, max)."""
    n = len(values)
    if n == 0:
        return EMPTY
    mean = math.fsum(values) / n
    m2 = math.fsum((v - mean) ** 2 for v in values)
    return (n, mean, m2, min(values), max(values))


def acc_merge(a, b):
    """Combine two accumulators without a second pass over the cells."""
    if a[0] == 0:
        return b
    if b[0] == 0:
        return a
    n = a[0] + b[0]
    delta = b[1] - a[1]
    mean = a[1] + delta * b[0] / n
    m2 = a[2] + b[2] + delta * delta * a[0] * b[0] / n
    return (n, mean, m2, min(a[3], b[3]), max(a[4], b[4]))


def summarize(acc):
    """min, max, mean and population standard deviation, or None when empty.

    GDAL divides by the cell count, not by count minus one.
    """
    if acc[0] == 0:
        return None
    return {"min": float(acc[3]), "max": float(acc[4]), "mean": float(acc[1]),
            "std": math.sqrt(acc[2] / acc[0])}


class Scan(object):
    """What one pass over the cells found."""

    def __init__(self, valid, nodata, nd_hits, nan_cells, total):
        self.valid = valid          # accumulator of the unmasked cells
        self.nodata = nodata        # declared NoData value or None
        self.nd_hits = nd_hits      # finite cells equal to the NoData value
        self.nan_cells = nan_cells  # NaN cells
        self.total = total

    @property
    def nodata_is_nan(self):
        return self.nodata is not None and self.nodata != self.nodata

    @property
    def nodata_cells(self):
        """Cells the declared NoData value masks. NaN is NoData when declared."""
        return self.nd_hits + (self.nan_cells if self.nodata_is_nan else 0)

    def unmasked(self):
        """Accumulator as if no NoData value were declared (NaN still out)."""
        if self.nd_hits == 0:
            return self.valid
        return acc_merge(self.valid, (self.nd_hits, self.nodata, 0.0,
                                      self.nodata, self.nodata))


def scan_cells(chunks, nodata, can_nan):
    """Fold an iterable of cell lists into a Scan. NaN never counts as a value."""
    valid = EMPTY
    hits = nans = total = 0
    mask = nodata is not None and nodata == nodata
    for chunk in chunks:
        total += len(chunk)
        if can_nan:
            kept = [v for v in chunk if v == v]
            nans += len(chunk) - len(kept)
            chunk = kept
        if mask:
            found = chunk.count(nodata)
            if found:
                hits += found
                chunk = [v for v in chunk if v != nodata]
        valid = acc_merge(valid, acc_of(chunk))
    return Scan(valid, nodata, hits, nans, total)


def close(a, b, tol):
    """Equal within a relative tolerance, with a floor of 1.0 on the scale."""
    return a == b or abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def valid_percent_text(valid, total):
    """The STATISTICS_VALID_PERCENT text GDAL writes for these counts."""
    if valid == 0:
        return "0"
    if valid == total:
        return "100"
    text = "%.4g" % (100.0 * valid / total)
    if text == "100":
        return "99.999"
    return text


def nodata_problem(nodata, int_range):
    """Why a declared NoData value can never equal a cell, or None."""
    if int_range is None or nodata is None:
        return None
    lo, hi = int_range
    if nodata != nodata or nodata in (float("inf"), float("-inf")):
        return "NoData %s is not a number an integer raster can hold" % nodata
    if nodata != math.floor(nodata) or nodata < lo or nodata > hi:
        return ("NoData %s is not an integer inside %d..%d, so no cell can "
                "equal it" % (_g(nodata), lo, hi))
    return None


def _g(x):
    return "%.10g" % x


def _mismatches(stored, summary, vpct, tol):
    """(name, stored, other) for each stored number that disagrees."""
    out = []
    for key, name in STAT_KEYS:
        if key in stored and not close(stored[key], summary[name], tol):
            out.append((name, stored[key], summary[name]))
    if K_VALID in stored and abs(stored[K_VALID] - vpct) > VALID_PCT_TOL:
        out.append(("valid_percent", stored[K_VALID], vpct))
    return out


def judge(scan, stored, approx, tol=DEFAULT_TOL, require_stats=False,
          problem=None):
    """Compare stored numbers with the cells. Returns (recomputed, findings)."""
    findings = []

    def add(code, severity, message):
        findings.append({"code": code, "severity": severity,
                         "message": message})

    if problem:
        add("NODATA_UNREPRESENTABLE", WARN, problem)
    n = scan.valid[0]
    if n == 0:
        add("ALL_NODATA", FAIL,
            "no valid cell: all %d cells are NoData or NaN, so no statistic "
            "exists" % scan.total)
        return None, findings
    rec = summarize(scan.valid)
    if not stored:
        add("NO_STATS", FAIL if require_stats else WARN,
            "no stored statistics: a legend or stretch has nothing to read")
    else:
        vp_rec = float(valid_percent_text(n, scan.total))
        bad = _mismatches(stored, rec, vp_rec, tol)
        if bad:
            nd = scan.nodata_cells
            vp_all = float(valid_percent_text(n + scan.nd_hits, scan.total))
            if scan.nd_hits and not _mismatches(
                    stored, summarize(scan.unmasked()), vp_all, tol):
                code = "STATS_COUNT_NODATA"
                text = ("stored statistics include the %d cells equal to the "
                        "declared NoData value %s: either the NoData value "
                        "changed after the statistics were stored, or it "
                        "masks real cells" % (nd, _g(scan.nodata)))
            else:
                code = "STALE_STATS"
                text = "stored statistics differ from the cells: " + "; ".join(
                    "%s stored %s, recomputed %s" % (name, _g(s), _g(r))
                    for name, s, r in bad)
            if approx:
                text += " (stored statistics are flagged approximate)"
            add(code, WARN if approx else FAIL, text)
    if scan.nd_hits and rec["min"] < scan.nodata < rec["max"]:
        add("NODATA_INSIDE_RANGE", WARN,
            "NoData value %s lies inside the data range %s..%s and masks %d "
            "cells (%.2f percent): confirm they are not real values"
            % (_g(scan.nodata), _g(rec["min"]), _g(rec["max"]), scan.nd_hits,
               100.0 * scan.nd_hits / scan.total))
    return rec, findings


def verdict_of(findings):
    severities = set(f["severity"] for f in findings)
    if FAIL in severities:
        return FAIL
    if WARN in severities:
        return WARN
    return OK


# ------------------------------------------------------------------- readers

class Raster(object):
    """A readable raster: size, NoData, a cell iterator and stored numbers."""

    def __init__(self, kind, width, height, nodata, can_nan, int_range, raw,
                 chunks):
        self.kind = kind
        self.width = width
        self.height = height
        self.nodata = nodata
        self.can_nan = can_nan
        self.int_range = int_range
        self.raw = raw          # stored metadata text by key
        self.chunks = chunks    # callable returning an iterator of cell lists


def _parse_float(text, what):
    try:
        return float(text)
    except ValueError:
        raise UnreadableRaster("%s %r is not a number" % (what, text[:40]))


def _asc_chunks(path, first, width, height):
    need = width * height
    seen = 0
    carry = b""
    with open(path, "rb") as handle:
        handle.seek(first)
        while True:
            block = handle.read(BLOCK_BYTES)
            data = carry + block
            tokens = data.split()
            if block and data[-1:].strip():
                carry = tokens.pop()
            else:
                carry = b""
            if tokens:
                try:
                    cells = [float(t) for t in tokens]
                except ValueError:
                    raise UnreadableRaster(
                        "a cell is not a number near cell %d" % (seen + 1))
                seen += len(cells)
                yield cells
            if not block:
                break
    if seen != need:
        raise UnreadableRaster("header promises %d cells, file holds %d"
                               % (need, seen))


def read_asc(path):
    header = {}
    first = 0
    with open(path, "rb") as handle:
        while True:
            line = handle.readline()
            parts = line.split()
            if not parts or parts[0].lower() not in ASC_KEYS:
                break
            if len(parts) != 2:
                raise UnreadableRaster("header line %r has no single value"
                                       % line[:40])
            header[parts[0].lower()] = parts[1]
            first = handle.tell()
    dims = []
    for key in (b"ncols", b"nrows"):
        if key not in header:
            raise UnreadableRaster("header has no %s" % key.decode("ascii"))
        try:
            dims.append(int(header[key]))
        except ValueError:
            raise UnreadableRaster("header %s is not an integer"
                                   % key.decode("ascii"))
    width, height = dims
    if width < 1 or height < 1:
        raise UnreadableRaster("header size %d x %d is empty" % (width, height))
    nodata = None
    if b"nodata_value" in header:
        nodata = _parse_float(header[b"nodata_value"].decode("latin-1"),
                              "header NODATA_value")
    return Raster("asc", width, height, nodata, True, None, {},
                  lambda: _asc_chunks(path, first, width, height))


# Sample format and bit depth to array typecode, kept only where the typecode
# has the width the file needs on this interpreter.
_DTYPES = dict(
    (key, code) for key, code in (
        ((1, 8), "B"), ((2, 8), "b"), ((1, 16), "H"), ((2, 16), "h"),
        ((1, 32), "I"), ((2, 32), "i"), ((3, 32), "f"), ((3, 64), "d"))
    if array.array(code).itemsize == key[1] // 8)

_TYPE_SIZE = {1: 1, 2: 1, 3: 2, 4: 4}
_TYPE_CODE = {1: "B", 3: "H", 4: "I"}
COMPRESSIONS = (1, 5, 8, 32946)


def _decode_tag(e, typ, cnt, data):
    if typ == 2:
        return data.split(b"\0")[0].decode("utf-8", "replace")
    return list(struct.unpack(e + str(cnt) + _TYPE_CODE[typ], data))


def _read_tags(handle, e, offset, size):
    handle.seek(offset)
    raw = handle.read(2)
    if len(raw) < 2:
        raise UnreadableRaster("the first IFD lies outside the file")
    count = struct.unpack(e + "H", raw)[0]
    body = handle.read(12 * count)
    if len(body) < 12 * count:
        raise UnreadableRaster("the IFD is truncated")
    tags = {}
    for i in range(count):
        tag, typ, cnt, field = struct.unpack(e + "HHI4s",
                                             body[12 * i:12 * i + 12])
        width = _TYPE_SIZE.get(typ)
        if width is None:
            continue
        nbytes = width * cnt
        if nbytes <= 4:
            data = field[:nbytes]
        else:
            where = struct.unpack(e + "I", field)[0]
            if where + nbytes > size:
                raise UnreadableRaster("tag %d points outside the file" % tag)
            handle.seek(where)
            data = handle.read(nbytes)
        tags[tag] = _decode_tag(e, typ, cnt, data)
    return tags


def lzw_decode(data):
    """TIFF LZW: MSB first, 9 to 12 bit codes, early change."""
    out = bytearray()
    table = None
    prev = None
    width = 9
    buf = bits = 0
    for byte in bytearray(data):
        buf = (buf << 8) | byte
        bits += 8
        while bits >= width:
            code = (buf >> (bits - width)) & ((1 << width) - 1)
            bits -= width
            buf &= (1 << bits) - 1
            if code == 257:
                return bytes(out)
            if code == 256:
                table = [bytes([i]) for i in range(256)]
                table += [b"", b""]
                width = 9
                prev = None
                continue
            if table is None:
                raise UnreadableRaster("LZW data does not start with a clear "
                                       "code")
            if prev is None:
                if code >= 256:
                    raise UnreadableRaster("LZW code %d before any entry"
                                           % code)
                entry = table[code]
            else:
                if code < len(table):
                    entry = table[code]
                elif code == len(table):
                    entry = prev + prev[:1]
                else:
                    raise UnreadableRaster("LZW code %d is not in the table"
                                           % code)
                table.append(prev + entry[:1])
            out += entry
            prev = entry
            size = len(table)
            width = 12 if size >= 2047 else 11 if size >= 1023 else (
                10 if size >= 511 else 9)
    return bytes(out)


def _cells(raw, code, swap):
    values = array.array(code)
    values.frombytes(raw)
    if swap:
        values.byteswap()
    return values.tolist()


def _undo_predictor(values, rowlen, bits, signed):
    mod = 1 << bits
    half = mod >> 1
    out = []
    for start in range(0, len(values), rowlen):
        for v in itertools.accumulate(values[start:start + rowlen]):
            v %= mod
            out.append(v - mod if signed and v >= half else v)
    return out


def _unit_values(handle, lay, off, cnt, ncells, rowlen, swap):
    """All cells of one compressed or predicted strip or tile."""
    handle.seek(off)
    raw = handle.read(cnt)
    if len(raw) < cnt:
        raise UnreadableRaster("the file is truncated inside a block")
    if lay["comp"] == 5:
        raw = lzw_decode(raw)
    elif lay["comp"] in (8, 32946):
        try:
            raw = zlib.decompress(raw)
        except zlib.error as exc:
            raise UnreadableRaster("deflate data is corrupt: %s" % exc)
    need = ncells * lay["isz"]
    if len(raw) < need:
        raise UnreadableRaster("a block decodes to %d bytes, %d needed"
                               % (len(raw), need))
    if lay["pred"] == 1:
        return _cells(raw[:need], lay["code"], swap)
    values = _cells(raw[:need], lay["code"].upper(), swap)
    return _undo_predictor(values, rowlen, lay["bits"], lay["code"].islower())


def _blocked(handle, lay, off, cnt, ncells, swap):
    """An uncompressed strip, read in bounded blocks."""
    need = ncells * lay["isz"]
    if cnt < need:
        raise UnreadableRaster("a strip holds %d bytes, %d needed"
                               % (cnt, need))
    handle.seek(off)
    step = max(1, BLOCK_BYTES // lay["isz"]) * lay["isz"]
    while need:
        take = min(step, need)
        raw = handle.read(take)
        if len(raw) < take:
            raise UnreadableRaster("the file is truncated inside a strip")
        yield _cells(raw, lay["code"], swap)
        need -= take


def _tiff_chunks(path, lay):
    swap = (lay["e"] == "<") != (sys.byteorder == "little")
    width, height = lay["width"], lay["height"]
    with open(path, "rb") as handle:
        if lay["tile"] is None:
            plain = lay["comp"] == 1 and lay["pred"] == 1
            for index, (off, cnt) in enumerate(lay["blocks"]):
                rows = min(lay["rps"], height - index * lay["rps"])
                if plain:
                    for chunk in _blocked(handle, lay, off, cnt, rows * width,
                                          swap):
                        yield chunk
                else:
                    yield _unit_values(handle, lay, off, cnt, rows * width,
                                       width, swap)
        else:
            tw, th = lay["tile"]
            across = -(-width // tw)
            for index, (off, cnt) in enumerate(lay["blocks"]):
                ty, tx = divmod(index, across)
                values = _unit_values(handle, lay, off, cnt, tw * th, tw, swap)
                vw = min(tw, width - tx * tw)
                vh = min(th, height - ty * th)
                if vw == tw and vh == th:
                    yield values
                else:
                    yield [v for r in range(vh)
                           for v in values[r * tw:r * tw + vw]]


def _first(tags, tag, default=None):
    value = tags.get(tag)
    if value is None:
        if default is None:
            raise UnreadableRaster("TIFF tag %d is missing" % tag)
        return default
    if not isinstance(value, list) or not value:
        raise UnreadableRaster("TIFF tag %d has the wrong type" % tag)
    return value[0]


def _ascii(tags, tag):
    value = tags.get(tag)
    if value is not None and not isinstance(value, str):
        raise UnreadableRaster("TIFF tag %d is not text" % tag)
    return value


def _as_float32(value):
    """Round to float32. Some Python versions raise past range, some give inf."""
    if abs(value) >= 3.4028235677973366e38:
        return math.copysign(float("inf"), value)
    return struct.unpack("f", struct.pack("f", value))[0]


def _tiff_layout(tags, e):
    lay = {"e": e}
    lay["width"] = _first(tags, 256)
    lay["height"] = _first(tags, 257)
    if lay["width"] < 1 or lay["height"] < 1:
        raise UnreadableRaster("TIFF size is empty")
    if _first(tags, 277, 1) != 1:
        raise UnreadableRaster("more than one band per pixel: nodatalie reads "
                               "single band rasters only")
    fmt = _first(tags, 339, 1)
    bits = _first(tags, 258, 1)
    code = _DTYPES.get((fmt, bits))
    if code is None:
        raise UnreadableRaster("sample format %d at %d bits is not supported"
                               % (fmt, bits))
    lay.update(code=code, bits=bits, isz=bits // 8, fmt=fmt)
    lay["comp"] = _first(tags, 259, 1)
    if lay["comp"] not in COMPRESSIONS:
        raise UnreadableRaster("TIFF compression %d is not supported"
                               % lay["comp"])
    lay["pred"] = _first(tags, 317, 1)
    if lay["pred"] not in (1, 2) or (
            lay["pred"] == 2 and (fmt == 3 or lay["comp"] == 1)):
        raise UnreadableRaster("TIFF predictor %d is not supported with this "
                               "sample format and compression" % lay["pred"])
    if 322 in tags:
        lay["tile"] = (_first(tags, 322), _first(tags, 323))
        offsets, counts = tags.get(324), tags.get(325)
        if lay["tile"][0] < 1 or lay["tile"][1] < 1:
            raise UnreadableRaster("TIFF tile size is empty")
        expect = (-(-lay["width"] // lay["tile"][0])
                  * -(-lay["height"] // lay["tile"][1]))
    else:
        lay["tile"] = None
        lay["rps"] = max(1, min(_first(tags, 278, lay["height"]),
                                lay["height"]))
        offsets, counts = tags.get(273), tags.get(279)
        expect = -(-lay["height"] // lay["rps"])
    if (not isinstance(offsets, list) or not isinstance(counts, list)
            or len(offsets) != expect or len(counts) != expect):
        raise UnreadableRaster("the block offsets do not match the image "
                               "layout (%d blocks expected)" % expect)
    lay["blocks"] = list(zip(offsets, counts))
    return lay


def _gdal_items(text):
    """STATISTICS_ items of band 0 from a GDAL_METADATA tag."""
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise UnreadableRaster("the GDAL metadata tag is not XML: %s" % exc)
    return dict((item.get("name"), (item.text or "").strip())
                for item in root.iter("Item") if item.get("sample") == "0")


def read_tiff(path):
    size = os.path.getsize(path)
    with open(path, "rb") as handle:
        head = handle.read(8)
        if head[:2] == b"II":
            e = "<"
        elif head[:2] == b"MM":
            e = ">"
        else:
            raise UnreadableRaster("not a TIFF file")
        if len(head) < 8:
            raise UnreadableRaster("the TIFF header is truncated")
        magic, ifd = struct.unpack(e + "HI", head[2:8])
        if magic == 43:
            raise UnreadableRaster("BigTIFF is not supported")
        if magic != 42:
            raise UnreadableRaster("not a TIFF file")
        tags = _read_tags(handle, e, ifd, size)
    lay = _tiff_layout(tags, e)
    nodata = None
    text = _ascii(tags, 42113)
    if text is not None:
        nodata = _parse_float(text.strip(), "GDAL_NODATA tag")
        if lay["fmt"] == 3 and lay["bits"] == 32:
            nodata = _as_float32(nodata)
    int_range = None
    if lay["fmt"] != 3:
        span = 1 << lay["bits"]
        int_range = (0, span - 1) if lay["fmt"] == 1 else (-span // 2,
                                                           span // 2 - 1)
    meta = _ascii(tags, 42112)
    raw = {} if meta is None else _gdal_items(meta)
    return Raster("tif", lay["width"], lay["height"], nodata, lay["fmt"] == 3,
                  int_range, raw, lambda: _tiff_chunks(path, lay))


def sidecar_path(path):
    return path + ".aux.xml"


def _band_one(root):
    for band in root.iter("PAMRasterBand"):
        if band.get("band") == "1":
            return band
    return None


def _mdi(band):
    out = {}
    for meta in band.findall("Metadata"):
        if meta.get("domain") is None:
            for item in meta.findall("MDI"):
                out[item.get("key")] = (item.text or "").strip()
    return out


def read_sidecar(path):
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise UnreadableRaster("the sidecar is not XML: %s" % exc)
    band = _band_one(root)
    return {} if band is None else _mdi(band)


def open_raster(path):
    """Open a raster and merge its stored numbers. The sidecar wins per key."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".asc":
        raster = read_asc(path)
    elif ext in (".tif", ".tiff"):
        raster = read_tiff(path)
    else:
        raise UnreadableRaster("%s files are not read" % (ext or "extensionless"))
    if os.path.exists(sidecar_path(path)):
        raster.raw.update(read_sidecar(sidecar_path(path)))
    return raster


def stored_numbers(raw):
    """(numbers by key, approximate flag) from the raw metadata text."""
    out = {}
    for key in [k for k, _n in STAT_KEYS] + [K_VALID]:
        if key in raw:
            out[key] = _parse_float(raw[key], "stored " + key)
    return out, K_APPROX in raw


# --------------------------------------------------------------- the checker

def check_file(path, tol=DEFAULT_TOL, require_stats=False):
    """Read one raster and return a report dict. Never raises on bad input."""
    try:
        raster = open_raster(path)
        stored, approx = stored_numbers(raster.raw)
        scan = scan_cells(raster.chunks(), raster.nodata, raster.can_nan)
    except (UnreadableRaster, OSError) as exc:
        return {"path": path, "verdict": UNREADABLE, "reason": str(exc),
                "findings": []}
    rec, findings = judge(scan, stored, approx, tol, require_stats,
                          nodata_problem(raster.nodata, raster.int_range))
    return {"path": path, "format": raster.kind, "verdict": verdict_of(findings),
            "cells": scan.total, "valid_cells": scan.valid[0],
            "nodata": raster.nodata, "nodata_cells": scan.nodata_cells,
            "nan_cells": scan.nan_cells,
            "stored": dict(stored),
            "recomputed": rec, "findings": findings}


def collect(paths):
    """Raster files under the given paths, plus the paths that do not exist."""
    files, missing = [], []
    for path in paths:
        if os.path.isdir(path):
            for dirpath, dirnames, names in os.walk(path):
                dirnames.sort()
                for name in sorted(names):
                    if name.lower().endswith(RASTER_EXT):
                        files.append(os.path.join(dirpath, name))
        elif os.path.isfile(path):
            files.append(path)
        else:
            missing.append(path)
    return files, missing


# ------------------------------------------------------------------- refresh

FIXABLE = ("STALE_STATS", "STATS_COUNT_NODATA", "NO_STATS")
BLOCKING = ("NODATA_INSIDE_RANGE", "NODATA_UNREPRESENTABLE")


def refusal(report):
    """Why --apply must not rewrite the statistics, or None when it may."""
    codes = set(f["code"] for f in report["findings"])
    if report["verdict"] == UNREADABLE:
        return "the raster could not be read"
    if not codes & set(FIXABLE):
        return "nothing to refresh"
    if codes & set(BLOCKING):
        return ("the NoData value needs a decision first; refreshing would "
                "hide it")
    return None


def _set_mdi(meta, key, value):
    for item in meta.findall("MDI"):
        if item.get("key") == key:
            item.text = value
            return
    ET.SubElement(meta, "MDI", {"key": key}).text = value


def refresh_sidecar(path, report):
    """Write the recomputed numbers into the PAM sidecar: (path, message).

    The old sidecar is copied to name.aux.xml.nodatalie.bak first, and the new
    one replaces it by rename, so a crash leaves one whole file.
    """
    side = sidecar_path(path)
    if os.path.exists(side):
        if os.path.exists(side + BAK_SUFFIX):
            return None, "a backup already exists: %s" % (side + BAK_SUFFIX)
        shutil.copy2(side, side + BAK_SUFFIX)
        tree = ET.parse(side)
    else:
        tree = ET.ElementTree(ET.Element("PAMDataset"))
    root = tree.getroot()
    band = _band_one(root)
    if band is None:
        band = ET.SubElement(root, "PAMRasterBand", {"band": "1"})
    meta = None
    for candidate in band.findall("Metadata"):
        if candidate.get("domain") is None:
            meta = candidate
    if meta is None:
        meta = ET.SubElement(band, "Metadata")
    rec = report["recomputed"]
    for key, name in STAT_KEYS:
        _set_mdi(meta, key, "%.14g" % rec[name])
    _set_mdi(meta, K_VALID,
             valid_percent_text(report["valid_cells"], report["cells"]))
    for item in meta.findall("MDI"):
        if item.get("key") == K_APPROX:
            meta.remove(item)
    ET.indent(tree, space="  ")
    handle, temp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(side)),
                                    suffix=".tmp")
    try:
        with os.fdopen(handle, "wb") as out:
            tree.write(out, encoding="utf-8")
        os.replace(temp, side)
    except OSError:
        with contextlib.suppress(OSError):
            os.remove(temp)
        raise
    return side, "refreshed"


def apply_one(path, report, tol, require_stats):
    """Refresh one raster's sidecar and re-check it. Returns the new report."""
    reason = refusal(report)
    if reason is None:
        try:
            side, reason = refresh_sidecar(path, report)
        except OSError as exc:
            side, reason = None, "could not write the sidecar: %s" % exc
        if side:
            again = check_file(path, tol, require_stats)
            again["apply"] = reason
            return again
    report["apply"] = "not refreshed: " + reason
    return report


# ----------------------------------------------------------------------- cli

_ORDER = {FAIL: 0, UNREADABLE: 1, WARN: 2, OK: 3}


def describe(report):
    """Plain text lines for one report."""
    if report["verdict"] == UNREADABLE:
        lines = ["%-10s %s  %s" % (UNREADABLE, report["path"],
                                   report["reason"])]
    else:
        nd = report["nodata"]
        lines = ["%-10s %s  cells=%d valid=%d nodata=%s nodata_cells=%d"
                 % (report["verdict"], report["path"], report["cells"],
                    report["valid_cells"], "none" if nd is None else _g(nd),
                    report["nodata_cells"])]
        for f in report["findings"]:
            lines.append("    %-22s %s" % (f["code"], f["message"]))
    if "apply" in report:
        lines.append("    apply: %s" % report["apply"])
    return lines


def exit_code(reports, missing, strict):
    verdicts = set(r["verdict"] for r in reports)
    if FAIL in verdicts or (strict and WARN in verdicts):
        return 1
    if UNREADABLE in verdicts or missing or not reports:
        return 2
    return 0


class _Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        sys.stderr.write("%s: error: %s\n" % (self.prog, message))
        sys.exit(64)


def _parse(argv):
    ap = _Parser(
        prog="nodatalie.py",
        allow_abbrev=False,
        description="Gate on stored versus recomputed raster statistics.",
        epilog="Nothing is written without --apply, and --apply writes only "
               "the PAM sidecar.",
    )
    ap.add_argument("paths", nargs="*",
                    help=".tif, .tiff or .asc files, or folders to search")
    ap.add_argument("--tol", type=float, default=DEFAULT_TOL,
                    help="relative tolerance for a stored number "
                         "(default %g)" % DEFAULT_TOL)
    ap.add_argument("--require-stats", dest="require_stats",
                    action="store_true",
                    help="a raster with no stored statistics fails, not warns")
    ap.add_argument("--strict", action="store_true",
                    help="a warning exits 1 as well")
    ap.add_argument("--json", action="store_true",
                    help="print one JSON document, worst raster first")
    ap.add_argument("--apply", action="store_true",
                    help="write the recomputed statistics into the sidecar of "
                         "each raster that fails only on stale or missing "
                         "statistics")
    ap.add_argument("--self-test", dest="self_test", action="store_true",
                    help="run the offline assertions and exit")
    return ap.parse_args(argv)


def main(argv=None):
    args = _parse(sys.argv[1:] if argv is None else argv)
    if args.self_test:
        return self_test()
    if not args.paths:
        print("error: give a raster or a folder. Use --self-test to verify "
              "the tool without data.", file=sys.stderr)
        return 64
    if args.tol < 0:
        print("error: --tol cannot be negative.", file=sys.stderr)
        return 64
    files, missing = collect(args.paths)
    for path in missing:
        print("error: no such path: %s" % path, file=sys.stderr)
    if not files and not missing:
        print("error: no .tif, .tiff or .asc file found.", file=sys.stderr)
    reports = [check_file(p, args.tol, args.require_stats) for p in files]
    if args.apply:
        reports = [apply_one(p, r, args.tol, args.require_stats)
                   for p, r in zip(files, reports)]
    reports.sort(key=lambda r: _ORDER[r["verdict"]])
    counts = dict((v, sum(1 for r in reports if r["verdict"] == v))
                  for v in _ORDER)
    summary = ("NODATALIE fail=%d warn=%d ok=%d unreadable=%d scanned=%d"
               % (counts[FAIL], counts[WARN], counts[OK], counts[UNREADABLE],
                  len(reports)))
    if args.json:
        print(json.dumps({"summary": summary, "files": reports}, indent=2,
                         sort_keys=True))
    else:
        for report in reports:
            for line in describe(report):
                print(line)
        print(summary)
    return exit_code(reports, missing, args.strict)


# ------------------------------------------------------------------ fixtures
# Writers for the self-test. They build synthetic rasters in pure Python, so
# the test needs no GDAL and no data.

def lzw_encode(data):
    """TIFF LZW encoder, used only to build test files."""
    out = bytearray()
    state = [0, 0]  # bit buffer, bit count

    def emit(code, width):
        state[0] = (state[0] << width) | code
        state[1] += width
        while state[1] >= 8:
            out.append((state[0] >> (state[1] - 8)) & 255)
            state[1] -= 8
        state[0] &= (1 << state[1]) - 1

    def fresh():
        return dict((bytes(bytearray([i])), i) for i in range(256)), 258, 9

    table, nxt, width = fresh()
    emit(256, width)
    word = b""
    for byte in bytearray(data):
        joined = word + bytes(bytearray([byte]))
        if joined in table:
            word = joined
            continue
        emit(table[word], width)
        table[joined] = nxt
        nxt += 1
        word = bytes(bytearray([byte]))
        if nxt >= 4093:
            emit(256, width)
            table, nxt, width = fresh()
        elif nxt > (1 << width) - 1:
            width += 1
    if word:
        emit(table[word], width)
        nxt += 1
        width += nxt > (1 << width) - 1
    emit(257, width)
    if state[1]:
        out.append((state[0] << (8 - state[1])) & 255)
    return bytes(out)


def _assemble(e, blobs, tags):
    body = b"".join(blobs)
    ifd_off = 8 + len(body)
    entries = sorted(tags.items())
    extra_off = ifd_off + 2 + 12 * len(entries) + 4
    extra = b""
    ent = b""
    for tag, (typ, vals) in entries:
        if typ == 2:
            data = vals.encode("ascii") + b"\0"
            cnt = len(data)
        else:
            cnt = len(vals)
            data = struct.pack(e + str(cnt) + {1: "B", 3: "H", 4: "I",
                                               7: "B"}[typ], *vals)
        if len(data) <= 4:
            field = data.ljust(4, b"\0")
        else:
            field = struct.pack(e + "I", extra_off + len(extra))
            extra += data
        ent += struct.pack(e + "HHI4s", tag, typ, cnt, field)
    head = (b"II" if e == "<" else b"MM") + struct.pack(e + "HI", 42, ifd_off)
    return (head + body + struct.pack(e + "H", len(entries)) + ent
            + struct.pack(e + "I", 0) + extra)


def _encode_unit(rows, code, e, comp, pred, bits):
    flat = []
    for row in rows:
        if pred == 2:
            mod = 1 << bits
            row = [row[0] % mod] + [(row[i] - row[i - 1]) % mod
                                    for i in range(1, len(row))]
        flat.extend(row)
    pack_code = code.upper() if pred == 2 else code
    raw = struct.pack(e + str(len(flat)) + pack_code, *flat)
    if comp == 5:
        return lzw_encode(raw)
    if comp == 8:
        return zlib.compress(raw)
    return raw


def make_tiff(path, rows, bits=16, fmt=2, endian="<", comp=1, pred=1,
              rps=None, tile=None, nodata=None, items=None, drop=(),
              extra=None):
    """Write a single band TIFF. items are GDAL_METADATA statistics."""
    h, w = len(rows), len(rows[0])
    code = {(1, 8): "B", (2, 8): "b", (1, 16): "H", (2, 16): "h",
            (1, 32): "I", (2, 32): "i", (3, 32): "f", (3, 64): "d"}[(fmt, bits)]
    tags = {256: (4, [w]), 257: (4, [h]), 258: (3, [bits]), 259: (3, [comp]),
            277: (3, [1]), 339: (3, [fmt])}
    if pred != 1:
        tags[317] = (3, [pred])
    units = []
    if tile is None:
        rps = rps or h
        for top in range(0, h, rps):
            units.append(rows[top:top + rps])
        tags[278] = (4, [rps])
    else:
        tw, th = tile
        for top in range(0, h, th):
            for left in range(0, w, tw):
                unit = []
                for y in range(top, top + th):
                    row = rows[y][left:left + tw] if y < h else []
                    unit.append(row + [0] * (tw - len(row)))
                units.append(unit)
        tags[322] = (4, [tw])
        tags[323] = (4, [th])
    blobs = [_encode_unit(u, code, endian, comp, pred, bits) for u in units]
    offsets, at = [], 8
    for blob in blobs:
        offsets.append(at)
        at += len(blob)
    counts = [len(b) for b in blobs]
    tags[324 if tile else 273] = (4, offsets)
    tags[325 if tile else 279] = (4, counts)
    if nodata is not None:
        tags[42113] = (2, nodata)
    if items:
        tags[42112] = (2, "<GDALMetadata>" + "".join(
            '<Item name="%s" sample="0">%s</Item>' % kv
            for kv in sorted(items.items())) + "</GDALMetadata>")
    for tag in drop:
        tags.pop(tag, None)
    tags.update(extra or {})
    with open(path, "wb") as handle:
        handle.write(_assemble(endian, blobs, tags))


def make_asc(path, rows, nodata="-9999", flat=False, header_extra=""):
    """Write an Esri ASCII grid. flat puts every cell on one line."""
    lines = ["ncols %d" % len(rows[0]), "nrows %d" % len(rows),
             "xllcorner 0", "yllcorner 0", "cellsize 1"]
    if nodata is not None:
        lines.append("NODATA_value %s" % nodata)
    body = [" ".join(str(v) for v in row) for row in rows]
    if flat:
        body = [" ".join(body)]
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + header_extra + "\n"
                     + "\n".join(body) + "\n")


def make_sidecar(path, items, approx=False):
    """Write name.aux.xml holding the given statistics."""
    side = sidecar_path(path)
    pairs = dict(items)
    if approx:
        pairs[K_APPROX] = "YES"
    body = "".join('<MDI key="%s">%s</MDI>' % kv for kv in sorted(pairs.items()))
    with open(side, "w") as handle:
        handle.write('<PAMDataset><PAMRasterBand band="1"><Metadata>%s'
                     "</Metadata></PAMRasterBand></PAMDataset>" % body)
    return side


def good_items(values, nodata=None, total=None):
    """Statistics a correct writer stores for the unmasked cells."""
    cells = [v for v in values if v == v and v != nodata]
    summary = summarize(acc_of(cells))
    items = {}
    for key, name in STAT_KEYS:
        items[key] = "%.14g" % summary[name]
    items[K_VALID] = valid_percent_text(len(cells), total or len(values))
    return items


# ------------------------------------------------------------------ self-test

class _Tally(object):
    """The check()/raises() harness behind --self-test."""

    def __init__(self):
        self.passed = 0
        self.failed = []

    def check(self, cond, label):
        if cond:
            self.passed += 1
            print("PASS  %s" % label)
        else:
            self.failed.append(label)
            print("FAIL  %s" % label)

    def raises(self, fn, label):
        try:
            fn()
        except UnreadableRaster:
            self.check(True, label)
        except Exception as exc:
            self.check(False, "%s (wrong exception %r)" % (label, exc))
        else:
            self.check(False, "%s (no error raised)" % label)

    def finish(self):
        print("-" * 68)
        total = self.passed + len(self.failed)
        if self.failed:
            print("%d assertions, %d failed" % (total, len(self.failed)))
            for label in self.failed:
                print("  FAILED: %s" % label)
            return 1
        print("%d assertions, 0 failed" % total)
        return 0


def self_test():
    """Assertions over the core, the readers and the CLI. No data needed."""
    tally = _Tally()
    check, raises = tally.check, tally.raises

    def approx(a, b, tol=1e-9):
        return abs(a - b) <= tol * max(1.0, abs(a), abs(b))

    def codes(report):
        return [f["code"] for f in report["findings"]]

    def found_codes(found):
        return [f["code"] for f in found]

    @contextlib.contextmanager
    def block_bytes(n):
        saved = globals()["BLOCK_BYTES"]
        globals()["BLOCK_BYTES"] = n
        try:
            yield
        finally:
            globals()["BLOCK_BYTES"] = saved

    def run(argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try:
                code = main(argv)
            except SystemExit as exc:
                code = exc.code
        return code, out.getvalue(), err.getvalue()

    def stats_of(scan):
        s = summarize(scan.valid)
        return {K_MIN: s["min"], K_MAX: s["max"], K_MEAN: s["mean"],
                K_STD: s["std"],
                K_VALID: float(valid_percent_text(scan.valid[0], scan.total))}

    def codes_stream(codes_list):
        """Pack 9 bit LZW codes, most significant bit first."""
        bits = "".join(format(c, "09b") for c in codes_list)
        bits += "0" * (-len(bits) % 8)
        return bytes(bytearray(int(bits[i:i + 8], 2)
                               for i in range(0, len(bits), 8)))

    print("nodatalie self-test: no GDAL, no arcpy, no network, no data")
    print("-" * 68)

    # ---- accumulators and the population standard deviation
    check(acc_of([]) == EMPTY, "no values give the empty accumulator")
    acc = acc_of([1, 2, 4, 6])
    check(acc[0] == 4 and acc[3] == 1 and acc[4] == 6, "count, min and max")
    check(approx(acc[1], 3.25) and approx(acc[2], 14.75), "mean and spread")
    check(summarize(EMPTY) is None, "an empty accumulator has no summary")
    summary = summarize(acc)
    check(approx(summary["std"], 1.9202864369672),
          "the deviation divides by the count, as GDAL does  <-- pinned defect")
    check(not approx(summary["std"], math.sqrt(14.75 / 3.0)),
          "the sample deviation is not what GDAL stores  <-- pinned defect")
    check(acc_merge(EMPTY, acc) is acc, "merging into empty returns the other")
    check(acc_merge(acc, EMPTY) is acc, "merging an empty changes nothing")
    whole = acc_of([5, 7, 9, 100, 3, 2])
    parts = acc_merge(acc_of([5, 7, 9]), acc_of([100, 3, 2]))
    check(parts[0] == 6 and approx(parts[1], whole[1])
          and approx(parts[2], whole[2]) and parts[3:] == (2, 100),
          "two halves merge to the whole")

    # ---- the scan: NoData, NaN and the unmasked view
    scan = scan_cells([[1, 2, -9999.0], [float("nan"), 4, 6]], -9999.0, True)
    check(scan.total == 6 and scan.valid[0] == 4, "six cells, four valid")
    check(scan.nd_hits == 1 and scan.nan_cells == 1, "one NoData, one NaN")
    check(scan.nodata_cells == 1, "NaN is not counted as the NoData value")
    check(scan.unmasked()[0] == 5 and scan.unmasked()[3] == -9999.0,
          "the unmasked view holds the NoData cell and not the NaN")
    nan_scan = scan_cells([[1, 2, float("nan")]], float("nan"), True)
    check(nan_scan.nodata_is_nan and nan_scan.nodata_cells == 1,
          "a NaN NoData value masks the NaN cells")
    check(nan_scan.valid[0] == 2 and nan_scan.nd_hits == 0,
          "a NaN NoData value finds no finite hit")
    check(nan_scan.unmasked() is nan_scan.valid,
          "with nothing to unmask the view is the valid accumulator")
    plain = scan_cells([[3, 4], [5]], None, False)
    check(plain.valid[0] == 3 and plain.nodata_cells == 0
          and not plain.nodata_is_nan, "no NoData value masks nothing")
    check(scan_cells([[1, 2], [-1, -1, -1]], -1, False).valid[0] == 2,
          "a NoData value is masked in every chunk")

    # ---- tolerance, percent, representability
    check(close(1.0, 1.0, 0.0), "equal numbers are close at zero tolerance")
    check(close(1000.0, 1000.0005, 1e-6), "a relative gap inside tolerance")
    check(not close(1000.0, 1000.01, 1e-6), "a relative gap outside tolerance")
    check(close(0.0, 1e-9, 1e-6), "the scale floor stops a tiny mean failing")
    check(close(float("inf"), float("inf"), 1e-6), "equal infinities match")
    check(not close(float("nan"), 1.0, 1e-6), "NaN matches nothing")
    check(valid_percent_text(0, 10) == "0", "no valid cell is 0 percent")
    check(valid_percent_text(10, 10) == "100", "all valid is 100 percent")
    check(valid_percent_text(2, 3) == "66.67", "two thirds is 66.67")
    check(valid_percent_text(99999, 100000) == "99.999",
          "never rounds a partly masked raster up to 100  <-- pinned defect")
    check(nodata_problem(-9999.0, None) is None, "floats take any NoData")
    check(nodata_problem(None, (0, 255)) is None, "no NoData has no problem")
    check(nodata_problem(0.0, (0, 255)) is None, "0 fits uint8")
    check(nodata_problem(255.0, (0, 255)) is None, "255 fits uint8")
    check(nodata_problem(-9999.0, (0, 255)) is not None,
          "-9999 cannot be a uint8 value")
    check(nodata_problem(256.0, (0, 255)) is not None, "256 is past uint8")
    check(nodata_problem(1.5, (-32768, 32767)) is not None,
          "a fraction cannot equal an integer cell")
    check(nodata_problem(float("nan"), (0, 255)) is not None,
          "NaN cannot equal an integer cell")
    check(nodata_problem(float("-inf"), (0, 255)) is not None,
          "infinity cannot equal an integer cell")

    # ---- judge: the gate itself
    cells = [10, 12, 14, 16, 11, 15, 17, 12, 14, 18]
    masked = cells[:1] + [-9999] + cells[1:5] + [-9999] + cells[5:]
    scan = scan_cells([masked], -9999, False)
    good = stats_of(scan)
    rec, found = judge(scan, good, False)
    check(found == [] and rec["min"] == 10.0, "matching statistics pass")
    check(verdict_of(found) == OK, "no finding is OK")
    rec, found = judge(scan, dict(good, STATISTICS_MINIMUM=-9999.0), False)
    check(found_codes(found) == ["STALE_STATS"]
          and found[0]["severity"] == FAIL, "a wrong minimum is stale")
    check("min stored -9999, recomputed 10" in found[0]["message"],
          "the finding names the number and both values")
    rec, found = judge(scan, {K_MIN: 10.0}, False)
    check(found == [], "a partial set that agrees passes")
    rec, found = judge(scan, {K_VALID: 50.0}, False)
    check(found_codes(found) == ["STALE_STATS"]
          and "valid_percent stored 50, recomputed 83.33" in found[0]["message"],
          "a valid percent that matches nothing is stale")
    rec, found = judge(scan, {K_VALID: 83.33}, False)
    check(found == [], "the stored valid percent is accepted at its rounding")
    unmasked = summarize(scan.unmasked())
    counted = {K_MIN: unmasked["min"], K_MAX: unmasked["max"],
               K_MEAN: unmasked["mean"], K_STD: unmasked["std"],
               K_VALID: 100.0}
    rec, found = judge(scan, counted, False)
    check(found_codes(found) == ["STATS_COUNT_NODATA"],
          "statistics that counted the NoData cells are named as such")
    check("2 cells equal to the declared NoData value -9999" in
          found[0]["message"], "that finding counts the cells")
    rec, found = judge(scan, {K_VALID: 100.0}, False)
    check(found_codes(found) == ["STATS_COUNT_NODATA"],
          "a valid percent of 100 over a masked raster counted the NoData")
    rec, found = judge(scan, counted, True)
    check(found[0]["severity"] == WARN
          and "flagged approximate" in found[0]["message"],
          "approximate statistics warn instead of failing")
    rec, found = judge(scan, dict(good, STATISTICS_STDDEV=5.0), True)
    check(found[0]["code"] == "STALE_STATS" and found[0]["severity"] == WARN,
          "a stale approximate set warns")
    sample_std = math.sqrt(acc_of([c for c in masked if c != -9999])[2] / 9.0)
    rec, found = judge(scan, dict(good, STATISTICS_STDDEV=sample_std), False)
    check(found_codes(found) == ["STALE_STATS"],
          "a sample deviation is not accepted as the stored one"
          "  <-- pinned defect")
    nudged = dict(good, STATISTICS_MEAN=good[K_MEAN] + 1e-4)
    rec, found = judge(scan, nudged, False, tol=0.0)
    check(found_codes(found) == ["STALE_STATS"],
          "a zero tolerance rejects a tiny difference")
    rec, found = judge(scan, nudged, False, tol=1e-3)
    check(found == [], "a loose tolerance accepts it")
    rec, found = judge(scan, {}, False)
    check(found_codes(found) == ["NO_STATS"] and found[0]["severity"] == WARN,
          "no stored statistics warn")
    rec, found = judge(scan, {}, False, require_stats=True)
    check(found[0]["severity"] == FAIL, "--require-stats makes that a failure")
    empty = scan_cells([[-9999, -9999]], -9999, False)
    rec, found = judge(empty, good, False, problem="NoData is odd")
    check(rec is None and found_codes(found) ==
          ["NODATA_UNREPRESENTABLE", "ALL_NODATA"],
          "an all NoData raster fails and keeps the NoData note")
    check(verdict_of(found) == FAIL, "FAIL outranks WARN")
    check(verdict_of(found[:1]) == WARN, "a lone warning is WARN")

    wide = scan_cells([[-5, 0, 0, 9, 40]], 0, False)
    rec, found = judge(wide, stats_of(wide), False)
    check(found_codes(found) == ["NODATA_INSIDE_RANGE"],
          "a NoData value inside the data range is flagged")
    check("masks 2 cells (40.00 percent)" in found[0]["message"],
          "that finding counts the masked cells")
    low = scan_cells([[0, 0, 3, 9]], 0, False)
    rec, found = judge(low, stats_of(low), False)
    check(found == [], "a NoData value below every real cell is fine")
    high = scan_cells([[-5, 3, 9, 9]], 9, False)
    rec, found = judge(high, stats_of(high), False)
    check(found == [], "a NoData value above every real cell is fine")
    unused = scan_cells([[1, 9]], 5, False)
    rec, found = judge(unused, stats_of(unused), False)
    check(found == [], "a NoData value no cell uses masks nothing")

    # ---- the LZW codec
    raw_stream = bytes(bytearray(range(256)) * 3)
    check(lzw_decode(lzw_encode(raw_stream)) == raw_stream,
          "the LZW codec round trips")
    check(lzw_decode(lzw_encode(b"a" * 50)) == b"a" * 50,
          "a code that is the entry being built is decoded")
    check(lzw_decode(lzw_encode(b"")) == b"", "LZW of nothing is nothing")
    check(lzw_decode(lzw_encode(b"abc")[:-1]) == b"abc",
          "a stream with no end code returns what it holds")
    raises(lambda: lzw_decode(b"\x7f\xff\xff"),
           "LZW data with no clear code is unreadable")
    raises(lambda: lzw_decode(codes_stream([256, 258])),
           "LZW with a table code before any entry is unreadable")
    raises(lambda: lzw_decode(codes_stream([256, 97, 300])),
           "LZW with a code outside the table is unreadable")

    # ---- fixtures on disk
    with tempfile.TemporaryDirectory() as tmp:
        def at(name):
            return os.path.join(tmp, name)

        def digest(path):
            with open(path, "rb") as handle:
                return handle.read()

        def put(name, text):
            with open(at(name), "w") as handle:
                handle.write(text)

        dem = [[10, 12, 14, 16], [11, -9999, 15, 17], [12, 14, -9999, 18]]
        flat_dem = [v for row in dem for v in row]
        items = good_items(flat_dem, -9999)

        # ---- the ASCII grid reader
        make_asc(at("ok.asc"), dem)
        make_sidecar(at("ok.asc"), items)
        rep = check_file(at("ok.asc"))
        check(rep["verdict"] == OK and rep["format"] == "asc",
              "a grid with correct statistics passes")
        check(rep["cells"] == 12 and rep["valid_cells"] == 10
              and rep["nodata_cells"] == 2, "its cells are counted")
        check(rep["nodata"] == -9999.0, "the header NoData value is read")
        make_asc(at("nostat.asc"), dem)
        check(codes(check_file(at("nostat.asc"))) == ["NO_STATS"],
              "a grid with no sidecar has no stored statistics")
        make_asc(at("nond.asc"), dem, nodata=None)
        rep = check_file(at("nond.asc"))
        check(rep["nodata"] is None and rep["valid_cells"] == 12,
              "a header with no NODATA_value declares none")
        make_asc(at("one.asc"), dem, flat=True)
        check(check_file(at("one.asc"))["cells"] == 12,
              "a grid on one line is read")
        with block_bytes(7):
            tiny = check_file(at("one.asc"))
        mean_valid = summarize(acc_of(cells))["mean"]
        check(tiny["cells"] == 12 and tiny["valid_cells"] == 10
              and approx(tiny["recomputed"]["mean"], mean_valid),
              "a token split across read blocks is rejoined"
              "  <-- pinned defect")
        put("nonl.asc", "ncols 2\nnrows 1\nnodata_value -9999\n3 4")
        check(check_file(at("nonl.asc"))["cells"] == 2,
              "a last cell with no newline is read  <-- pinned defect")
        put("upper.asc", "NCOLS 2\nNROWS 1\nXLLCENTER 0\nYLLCENTER 0\n"
            "DX 1\nDY 1\nNODATA_VALUE -1\n3 -1\n")
        rep = check_file(at("upper.asc"))
        check(rep["valid_cells"] == 1 and rep["nodata"] == -1.0,
              "upper case keys and centre corners are read")
        put("nan.asc", "ncols 3\nnrows 1\nnodata_value -9999\n1 NaN 3\n")
        rep = check_file(at("nan.asc"))
        check(rep["nan_cells"] == 1 and rep["valid_cells"] == 2,
              "a NaN cell in a grid is excluded and counted")
        put("short.asc", "ncols 3\nnrows 2\nnodata_value -9999\n1 2 3 4\n")
        rep = check_file(at("short.asc"))
        check(rep["verdict"] == UNREADABLE
              and "promises 6 cells, file holds 4" in rep["reason"],
              "a truncated grid is unreadable, never clean")
        put("long.asc", "ncols 1\nnrows 1\n1 2\n")
        check("file holds 2" in check_file(at("long.asc"))["reason"],
              "a grid with extra cells is unreadable")
        put("word.asc", "ncols 2\nnrows 1\n1 abc\n")
        check("not a number" in check_file(at("word.asc"))["reason"],
              "a cell that is not a number is unreadable")
        put("nocols.asc", "nrows 1\n1 2\n")
        check("no ncols" in check_file(at("nocols.asc"))["reason"],
              "a header with no ncols is unreadable")
        put("badcols.asc", "ncols x\nnrows 1\n1 2\n")
        check("not an integer" in check_file(at("badcols.asc"))["reason"],
              "an ncols that is not an integer is unreadable")
        put("zero.asc", "ncols 0\nnrows 1\n")
        check("empty" in check_file(at("zero.asc"))["reason"],
              "a zero width grid is unreadable")
        put("novalue.asc", "ncols\nnrows 1\n1 2\n")
        check("no single value" in check_file(at("novalue.asc"))["reason"],
              "a header key with no value is unreadable")
        put("badnd.asc", "ncols 1\nnrows 1\nnodata_value x\n1\n")
        check("NODATA_value" in check_file(at("badnd.asc"))["reason"],
              "a NoData value that is not a number is unreadable")
        put("empty.asc", "")
        check("no ncols" in check_file(at("empty.asc"))["reason"],
              "an empty file is unreadable")
        put("blank.asc", "ncols 2\nnrows 1\n\n1 2\n")
        check(check_file(at("blank.asc"))["cells"] == 2,
              "a blank line before the cells is skipped")

        # ---- the stored numbers: sidecar, tag, precedence
        make_asc(at("stale.asc"), dem)
        make_sidecar(at("stale.asc"), dict(items, STATISTICS_MINIMUM="-9999"))
        rep = check_file(at("stale.asc"))
        check(rep["verdict"] == FAIL and codes(rep) == ["STALE_STATS"],
              "a sidecar minimum that counted a sentinel fails")
        make_asc(at("count.asc"), dem)
        make_sidecar(at("count.asc"), good_items(flat_dem))
        rep = check_file(at("count.asc"))
        check(codes(rep) == ["STATS_COUNT_NODATA"] and rep["verdict"] == FAIL,
              "NoData added after the statistics were stored is named"
              "  <-- pinned defect")
        make_asc(at("approx.asc"), dem)
        make_sidecar(at("approx.asc"), dict(items, STATISTICS_MEAN="1"),
                     approx=True)
        check(check_file(at("approx.asc"))["verdict"] == WARN,
              "approximate stale statistics warn")
        make_asc(at("other.asc"), dem)
        put("other.asc.aux.xml",
            '<PAMDataset><PAMRasterBand band="2"><Metadata>'
            '<MDI key="STATISTICS_MINIMUM">1</MDI></Metadata>'
            '</PAMRasterBand></PAMDataset>')
        check(codes(check_file(at("other.asc"))) == ["NO_STATS"],
              "statistics stored for band 2 are not band 1's")
        make_asc(at("dom.asc"), dem)
        put("dom.asc.aux.xml",
            '<PAMDataset><PAMRasterBand band="1">'
            '<Metadata domain="IMAGE_STRUCTURE">'
            '<MDI key="STATISTICS_MINIMUM">1</MDI></Metadata>'
            '</PAMRasterBand></PAMDataset>')
        check(codes(check_file(at("dom.asc"))) == ["NO_STATS"],
              "a metadata block with a domain is not the statistics block")
        make_asc(at("nostat2.asc"), dem)
        put("nostat2.asc.aux.xml", "<PAMDataset/>")
        check(codes(check_file(at("nostat2.asc"))) == ["NO_STATS"],
              "a sidecar with no band holds no statistics")
        make_asc(at("bad.asc"), dem)
        put("bad.asc.aux.xml", "<PAMDataset>")
        check("not XML" in check_file(at("bad.asc"))["reason"],
              "a sidecar that is not XML is unreadable")
        make_asc(at("word2.asc"), dem)
        make_sidecar(at("word2.asc"), {K_MIN: "abc"})
        check("stored STATISTICS_MINIMUM" in
              check_file(at("word2.asc"))["reason"],
              "a stored number that is not a number is unreadable")
        make_asc(at("nan2.asc"), dem)
        make_sidecar(at("nan2.asc"), {K_MEAN: "nan"})
        check(codes(check_file(at("nan2.asc"))) == ["STALE_STATS"],
              "a stored NaN never matches")

        # ---- the TIFF reader, every sample type, both byte orders
        grid = [[(x * 7 + y * 3) % 50 - 10 for x in range(37)]
                for y in range(29)]
        want = summarize(acc_of([v for row in grid for v in row]))
        for fmt, bits in ((1, 8), (2, 8), (1, 16), (2, 16), (1, 32), (2, 32),
                          (3, 32), (3, 64)):
            rows = grid if fmt != 1 else [[abs(v) for v in r] for r in grid]
            exp = summarize(acc_of([v for r in rows for v in r]))
            for endian in ("<", ">"):
                name = at("t_%d_%d_%s.tif" % (fmt, bits, "le" if endian == "<" else "be"))
                make_tiff(name, rows, bits=bits, fmt=fmt, endian=endian)
                rep = check_file(name)
                check(rep["format"] == "tif" and rep["cells"] == 29 * 37
                      and all(approx(rep["recomputed"][k], exp[k])
                              for k in exp),
                      "sample format %d at %d bits, byte order %s, reads "
                      "exactly" % (fmt, bits, endian))
        layouts = (("strips of 1 row", dict(rps=1)),
                   ("strips of 4 rows", dict(rps=4)),
                   ("one strip", dict()),
                   ("tiles of 16", dict(tile=(16, 16))),
                   ("tiles that fit", dict(tile=(37, 29))))
        for comp, pred in ((1, 1), (5, 1), (5, 2), (8, 1), (8, 2)):
            for label, kw in layouts:
                name = at("l_%d_%d.tif" % (comp, pred))
                make_tiff(name, grid, bits=16, fmt=2, comp=comp, pred=pred,
                          **kw)
                rep = check_file(name)
                check(rep["cells"] == 29 * 37
                      and approx(rep["recomputed"]["mean"], want["mean"])
                      and approx(rep["recomputed"]["std"], want["std"]),
                      "compression %d predictor %d, %s" % (comp, pred, label))
        make_tiff(at("be_pred.tif"), grid, bits=32, fmt=2, endian=">", comp=8,
                  pred=2, tile=(16, 16))
        check(approx(check_file(at("be_pred.tif"))["recomputed"]["mean"],
                     want["mean"]),
              "big endian 32 bit with predictor 2 and tiles reads exactly")
        make_tiff(at("sneg.tif"), [[-3, -2, -1, 0, 1]], bits=8, fmt=2, comp=5,
                  pred=2)
        check(check_file(at("sneg.tif"))["recomputed"]["min"] == -3.0,
              "a signed 8 bit predictor wraps to the right sign")
        make_tiff(at("blocks.tif"), grid, bits=16, fmt=2, rps=29)
        with block_bytes(10):
            blocks = check_file(at("blocks.tif"))
        check(blocks["cells"] == 29 * 37
              and approx(blocks["recomputed"]["mean"], want["mean"]),
              "a strip larger than the read block is read in pieces")
        make_tiff(at("blocks1.tif"), [[1, 2], [3, 4]], bits=32, fmt=2)
        with block_bytes(1):
            small = check_file(at("blocks1.tif"))
        check(small["cells"] == 4,
              "a read block smaller than one cell still advances"
              "  <-- pinned defect")
        lzw_rows = [[(x * 31 + y * 17) % 251 for x in range(120)]
                    for y in range(90)]
        make_tiff(at("lzwbig.tif"), lzw_rows, bits=8, fmt=1, comp=5)
        want_lzw = summarize(acc_of([v for r in lzw_rows for v in r]))
        check(approx(check_file(at("lzwbig.tif"))["recomputed"]["mean"],
                     want_lzw["mean"]),
              "LZW across the 9 to 12 bit code widths  <-- pinned defect")
        seed = [12345]

        def rnd():
            seed[0] = (seed[0] * 1103515245 + 12345) & 0x7FFFFFFF
            return (seed[0] >> 16) & 255
        noise = [[rnd() for _x in range(90)] for _y in range(90)]
        make_tiff(at("lzwnoise.tif"), noise, bits=8, fmt=1, comp=5)
        want_noise = summarize(acc_of([v for r in noise for v in r]))
        check(approx(check_file(at("lzwnoise.tif"))["recomputed"]["mean"],
                     want_noise["mean"]),
              "LZW on noise that fills the table and forces a reset"
              "  <-- pinned defect")

        # ---- NoData tag and stored tag statistics
        rows = [[3, 4, 0, 6], [0, 8, 9, 10]]
        make_tiff(at("nd.tif"), rows, bits=16, fmt=1, nodata="0")
        rep = check_file(at("nd.tif"))
        check(rep["nodata"] == 0.0 and rep["nodata_cells"] == 2
              and rep["valid_cells"] == 6, "the GDAL_NODATA tag is read")
        make_tiff(at("nd_long.tif"), rows, bits=16, fmt=2,
                  nodata="-32768.0000000000000")
        check(check_file(at("nd_long.tif"))["nodata"] == -32768.0,
              "a NoData text too long to sit inline is read")
        fl = [[1.5, 2.5, -3.4028234663852886e38], [4.5, 5.5, 6.5]]
        make_tiff(at("f32nd.tif"), fl, bits=32, fmt=3,
                  nodata="-3.40282346638529e+38")
        rep = check_file(at("f32nd.tif"))
        check(rep["nodata_cells"] == 1 and rep["valid_cells"] == 5,
              "a 15 digit float32 NoData text still matches the cells"
              "  <-- pinned defect")
        make_tiff(at("f32big.tif"), fl, bits=32, fmt=3, nodata="1e300")
        check(check_file(at("f32big.tif"))["nodata"] == float("inf")
              and check_file(at("f32big.tif"))["nodata_cells"] == 0,
              "a NoData past float32 range becomes infinity and matches "
              "nothing")
        make_tiff(at("f64nd.tif"), fl, bits=64, fmt=3,
                  nodata="-3.4028234663852886e38")
        check(check_file(at("f64nd.tif"))["nodata_cells"] == 1,
              "a float64 NoData is compared as it is")
        fnan = [[1.0, float("nan"), 3.0], [4.0, float("nan"), 6.0]]
        make_tiff(at("nannd.tif"), fnan, bits=32, fmt=3, nodata="nan")
        rep = check_file(at("nannd.tif"))
        check(rep["nodata_cells"] == 2 and rep["valid_cells"] == 4,
              "a NaN NoData tag masks the NaN cells")
        make_tiff(at("nan_only.tif"), fnan, bits=32, fmt=3)
        rep = check_file(at("nan_only.tif"))
        check(rep["nan_cells"] == 2 and rep["valid_cells"] == 4
              and rep["nodata"] is None,
              "NaN cells are excluded even with no NoData value")
        make_tiff(at("u8bad.tif"), [[1, 2], [3, 4]], bits=8, fmt=1,
                  nodata="-9999")
        check("NODATA_UNREPRESENTABLE" in codes(check_file(at("u8bad.tif"))),
              "a NoData value an 8 bit raster cannot hold is flagged")
        tag_good = good_items([3, 4, 0, 6, 0, 8, 9, 10], 0)
        make_tiff(at("tagstat.tif"), rows, bits=16, fmt=1, nodata="0",
                  items=tag_good)
        check(check_file(at("tagstat.tif"))["verdict"] == OK,
              "statistics stored in the file's own tag are read")
        make_tiff(at("tagstale.tif"), rows, bits=16, fmt=1, nodata="0",
                  items=dict(tag_good, STATISTICS_MAXIMUM="99"))
        check(check_file(at("tagstale.tif"))["verdict"] == FAIL,
              "stale statistics in the tag fail")
        make_sidecar(at("tagstale.tif"), tag_good)
        check(check_file(at("tagstale.tif"))["verdict"] == OK,
              "a sidecar overrides the tag, as GDAL does  <-- pinned defect")
        make_tiff(at("tagok.tif"), rows, bits=16, fmt=1, nodata="0",
                  items=tag_good)
        make_sidecar(at("tagok.tif"), {K_MAX: "99"})
        check(codes(check_file(at("tagok.tif"))) == ["STALE_STATS"],
              "a stale sidecar key beats a correct tag")
        make_tiff(at("band2.tif"), rows, bits=16, fmt=1,
                  extra={42112: (2, '<GDALMetadata><Item name="STATISTICS_MAX'
                                 'IMUM" sample="1">99</Item></GDALMetadata>')})
        check(codes(check_file(at("band2.tif"))) == ["NO_STATS"],
              "tag statistics for another band are not read")
        make_tiff(at("badxml.tif"), rows, bits=16, fmt=1,
                  extra={42112: (2, "<GDALMetadata")})
        check("not XML" in check_file(at("badxml.tif"))["reason"],
              "a metadata tag that is not XML is unreadable")
        make_tiff(at("emptyitem.tif"), rows, bits=16, fmt=1,
                  extra={42112: (2, '<GDALMetadata><Item name="STATISTICS_MI'
                                 'NIMUM" sample="0"/></GDALMetadata>')})
        check("not a number" in check_file(at("emptyitem.tif"))["reason"],
              "an empty stored value is unreadable")
        make_tiff(at("utype.tif"), rows, bits=16, fmt=1,
                  extra={50000: (7, [1, 2, 3, 4, 5, 6, 7])})
        check(check_file(at("utype.tif"))["cells"] == 8,
              "a tag of a type the reader does not need is skipped")
        make_tiff(at("ndtext.tif"), rows, bits=16, fmt=1,
                  extra={42113: (3, [5])})
        check("is not text" in check_file(at("ndtext.tif"))["reason"],
              "a NoData tag that is not text is unreadable")
        make_tiff(at("ndword.tif"), rows, bits=16, fmt=1, nodata="none")
        check("GDAL_NODATA" in check_file(at("ndword.tif"))["reason"],
              "a NoData tag that is not a number is unreadable")

        # ---- unreadable TIFFs, each with its reason
        def broken(name, mutate, needle, label, **kw):
            path = at(name)
            make_tiff(path, kw.pop("rows", rows),
                      **dict(dict(bits=16, fmt=1), **kw))
            data = bytearray(digest(path))
            mutate(data)
            with open(path, "wb") as handle:
                handle.write(bytes(data))
            check(needle in check_file(path).get("reason", ""), label)

        def fat_ifd(data):
            where = struct.unpack("<I", bytes(data[4:8]))[0]
            data[where:where + 2] = b"\xff\xff"

        def noop(data):
            return None

        broken("e1.tif", lambda d: d.__setitem__(slice(0, 2), b"XX"),
               "not a TIFF", "a file that is not a TIFF is unreadable")
        broken("e2.tif", lambda d: d.__setitem__(slice(2, 4), b"\x2b\x00"),
               "BigTIFF", "BigTIFF is named")
        broken("e3.tif", lambda d: d.__setitem__(slice(2, 4), b"\x07\x00"),
               "not a TIFF", "a wrong magic number is unreadable")
        broken("e4.tif", lambda d: d.__delitem__(slice(3, None)),
               "truncated", "a header cut short is unreadable")
        broken("e5.tif",
               lambda d: d.__setitem__(slice(4, 8), b"\xff\xff\xff\x7f"),
               "outside the file", "an IFD offset past the end is unreadable")
        broken("e6.tif", fat_ifd, "truncated",
               "an IFD that promises more entries than it holds is unreadable")
        broken("e7.tif", lambda d: d.__delitem__(slice(len(d) - 3, None)),
               "outside the file",
               "a tag whose data lies past the end of the file is unreadable",
               rows=[[1, 2], [3, 4], [5, 6]], rps=1)
        put("zerobytes.tif", "")
        check("not a TIFF" in check_file(at("zerobytes.tif"))["reason"],
              "an empty TIFF is unreadable")
        broken("e8.tif", noop, "missing",
               "a missing width tag is unreadable", drop=(256,))
        broken("e9.tif", noop, "wrong type",
               "a width stored as text is unreadable", extra={256: (2, "4")})
        broken("e10.tif", noop, "wrong type",
               "an empty width is unreadable", extra={256: (3, [])})
        broken("e11.tif", noop, "size is empty",
               "a zero height is unreadable", extra={257: (4, [0])})
        broken("e12.tif", noop, "single band",
               "three samples per pixel is unreadable", extra={277: (3, [3])})
        broken("e13.tif", noop, "not supported",
               "12 bit samples are unreadable", extra={258: (3, [12])})
        broken("e14.tif", noop, "compression 7",
               "an unsupported compression is named", extra={259: (3, [7])})
        broken("e15.tif", noop, "predictor 3",
               "predictor 3 is named", extra={317: (3, [3])}, comp=8)
        broken("e16.tif", noop, "predictor 2",
               "predictor 2 on floats is refused", comp=8,
               extra={317: (3, [2])}, fmt=3, bits=32)
        broken("e17.tif", noop, "predictor 2",
               "predictor 2 with no compression is refused  <-- pinned defect",
               extra={317: (3, [2])})
        broken("e18.tif", noop, "do not match",
               "a strip count that disagrees with the height is unreadable",
               extra={273: (4, [8, 8, 8])})
        broken("e19.tif", noop, "do not match",
               "tiles with no offsets are unreadable", tile=(16, 16),
               drop=(324,))
        broken("e20.tif", noop, "tile size is empty",
               "a zero tile width is unreadable", tile=(16, 16),
               extra={322: (4, [0])})
        broken("e21.tif", noop, "missing",
               "a tile width with no tile length is unreadable",
               tile=(16, 16), drop=(323,))
        broken("e22.tif", noop, "holds 10 bytes",
               "a byte count smaller than the strip is unreadable",
               extra={279: (4, [10])})
        broken("e23.tif", noop, "decodes to",
               "a block that decodes short is unreadable", comp=8,
               rows=[[1, 2, 3, 4]], extra={257: (4, [9])}, drop=(278,))

        def near_end(name, **kw):
            """A file whose only block starts ten bytes before its end."""
            probe = at(name)
            make_tiff(probe, grid, **kw)
            size = os.path.getsize(probe)
            make_tiff(probe, grid, extra={273: (4, [size - 10])}, **kw)
            return probe

        check("truncated inside a strip" in check_file(
            near_end("trunc_u.tif", bits=16, fmt=2))["reason"],
            "a strip with its bytes missing is unreadable, never clean")
        check("truncated inside a block" in check_file(
            near_end("trunc_c.tif", bits=16, fmt=2, comp=8))["reason"],
            "a compressed block with its bytes missing is unreadable")
        make_tiff(at("zbad.tif"), grid, bits=16, fmt=2, comp=8)
        zdata = bytearray(digest(at("zbad.tif")))
        for i in range(12, 40):
            zdata[i] ^= 0x5A
        with open(at("zbad.tif"), "wb") as handle:
            handle.write(bytes(zdata))
        check("deflate data is corrupt" in check_file(at("zbad.tif"))["reason"],
              "corrupt deflate data is unreadable")

        # ---- other files
        put("x.png", "x")
        check(".png files are not read" in check_file(at("x.png"))["reason"],
              "an unsupported extension is unreadable")
        put("noext", "x")
        check("extensionless" in check_file(at("noext"))["reason"],
              "a file with no extension is unreadable")
        check(check_file(at("missing.tif"))["verdict"] == UNREADABLE,
              "a file that does not exist is unreadable, never clean")

        # ---- walking a tree
        tree = at("tree")
        os.makedirs(os.path.join(tree, "b"))
        os.makedirs(os.path.join(tree, "a"))
        for rel in ("b/two.TIF", "a/one.asc", "a/skip.txt", "top.tiff",
                    "a/one.asc.aux.xml"):
            with open(os.path.join(tree, rel), "wb") as handle:
                handle.write(b"x")
        files, missing = collect([tree, at("nope")])
        names = [os.path.relpath(f, tree).replace(os.sep, "/") for f in files]
        check(names == ["top.tiff", "a/one.asc", "b/two.TIF"],
              "a tree yields rasters only, sorted within each folder")
        check(missing == [at("nope")], "a missing path is reported")
        files, missing = collect([at("ok.asc")])
        check(files == [at("ok.asc")] and missing == [],
              "a file argument is taken as it is")

        # ---- the command line
        good_dir = at("cli_ok")
        os.makedirs(good_dir)
        make_asc(os.path.join(good_dir, "g.asc"), dem)
        make_sidecar(os.path.join(good_dir, "g.asc"), items)
        code, out, err = run([good_dir])
        check(code == 0 and "NODATALIE fail=0 warn=0 ok=1" in out,
              "a clean tree exits 0")
        check(out.splitlines()[0].startswith("OK         "),
              "the report line starts with the verdict")
        bad_dir = at("cli_bad")
        os.makedirs(bad_dir)
        for name in ("a_stale.asc", "b_nostat.asc", "c_ok.asc"):
            make_asc(os.path.join(bad_dir, name), dem)
        make_sidecar(os.path.join(bad_dir, "a_stale.asc"),
                     dict(items, STATISTICS_MEAN="1"))
        make_sidecar(os.path.join(bad_dir, "c_ok.asc"), items)
        with open(os.path.join(bad_dir, "d_bad.asc"), "w") as handle:
            handle.write("ncols 9\n")
        code, out, err = run([bad_dir])
        lines = out.splitlines()
        check(code == 1, "a stale raster exits 1")
        check([ln.split()[0] for ln in lines if not ln.startswith(" ")][:4] ==
              ["FAIL", "UNREADABLE", "WARN", "OK"],
              "the worst verdict is listed first")
        check("fail=1 warn=1 ok=1 unreadable=1 scanned=4" in lines[-1],
              "the summary counts every verdict")
        check(any(ln.startswith("    STALE_STATS") for ln in lines),
              "the finding is printed under its raster")
        code, out, err = run([bad_dir, "--json"])
        doc = json.loads(out)
        check(doc["files"][0]["verdict"] == FAIL
              and doc["files"][-1]["verdict"] == OK, "--json lists worst first")
        check(doc["files"][0]["findings"][0]["code"] == "STALE_STATS"
              and doc["summary"].startswith("NODATALIE"),
              "--json carries the findings and the summary")
        check(run([good_dir, "--strict"])[0] == 0,
              "--strict leaves a clean tree at 0")
        warn_dir = at("cli_warn")
        os.makedirs(warn_dir)
        make_asc(os.path.join(warn_dir, "w.asc"), dem)
        check(run([warn_dir])[0] == 0, "a warning alone exits 0")
        check(run([warn_dir, "--strict"])[0] == 1, "--strict makes it 1")
        check(run([warn_dir, "--require-stats"])[0] == 1,
              "--require-stats makes a missing set fail")
        check(run([os.path.join(bad_dir, "d_bad.asc")])[0] == 2,
              "an unreadable raster alone exits 2")
        code, out, err = run([good_dir, at("nope")])
        check(code == 2 and "no such path" in err,
              "a path that does not exist exits 2, never 0")
        check(run([bad_dir, at("nope")])[0] == 1,
              "a failure outranks a missing path")
        empty_dir = at("cli_empty")
        os.makedirs(empty_dir)
        code, out, err = run([empty_dir])
        check(code == 2 and "no .tif" in err,
              "a folder with no raster exits 2, never 0")
        code, out, err = run([])
        check(code == 64 and "--self-test" in err,
              "no path is a usage error that names --self-test")
        check(run([good_dir, "--tol", "-1"])[0] == 64,
              "a negative tolerance is a usage error")
        check(run([good_dir, "--tol", "x"])[0] == 64,
              "a tolerance that is not a number is a usage error")
        check(run([good_dir, "--tol", "0.5"])[0] == 0,
              "a loose tolerance is accepted")
        code, out, err = run([good_dir, "--ap"])
        check(code == 64 and "unrecognized" in err,
              "a unique prefix --ap of --apply is REFUSED by the parser"
              "  <-- pinned defect")
        check(run([good_dir, "--app"])[0] == 64,
              "a longer prefix --app of --apply is REFUSED  <-- pinned defect")
        check(run([good_dir, "--self-t"])[0] == 64,
              "a unique prefix of --self-test is REFUSED  <-- pinned defect")
        check(run([good_dir, "--js"])[0] == 64,
              "a unique prefix of --json is REFUSED  <-- pinned defect")

        # ---- --apply, the only writer
        wdir = at("apply")
        os.makedirs(wdir)

        def wat(name):
            return os.path.join(wdir, name)

        wrows = [[10, 12, 14, 16], [11, -9999, 15, 17], [12, 14, -9999, 18]]
        for name in ("none.tif", "stale.tif", "clean.tif", "nobandsidecar.tif",
                     "nobandfile.tif"):
            make_tiff(wat(name), wrows, bits=16, fmt=2, nodata="-9999")
        old_side = wat("stale.tif.aux.xml")
        with open(old_side, "w") as handle:
            handle.write(
                '<PAMDataset><PAMRasterBand band="1"><Histograms>'
                '<HistItem><Min>0</Min></HistItem></Histograms>'
                '<Metadata domain="X"><MDI key="K">v</MDI></Metadata>'
                '<Metadata><MDI key="STATISTICS_MAXIMUM">999</MDI>'
                '<MDI key="STATISTICS_APPROXIMATE">YES</MDI>'
                '<MDI key="EXTRA_KEY">keep me</MDI></Metadata>'
                '</PAMRasterBand></PAMDataset>')
        old_bytes = digest(old_side)
        raster_before = digest(wat("stale.tif"))
        make_tiff(wat("inside.tif"), [[-5, 0, 0, 9, 40]], bits=16, fmt=2,
                  nodata="0")
        make_tiff(wat("alln.tif"), [[0, 0]], bits=16, fmt=2, nodata="0")
        make_tiff(wat("u8bad.tif"), [[1, 2]], bits=8, fmt=1, nodata="-9999")
        make_sidecar(wat("clean.tif"), items)
        with open(wat("nobandsidecar.tif.aux.xml"), "w") as handle:
            handle.write('<PAMDataset><PAMRasterBand band="1"><Metadata '
                         'domain="Q"><MDI key="a">b</MDI></Metadata>'
                         '</PAMRasterBand></PAMDataset>')
        with open(wat("nobandfile.tif.aux.xml"), "w") as handle:
            handle.write("<PAMDataset><Other/></PAMDataset>")
        with open(wat("broken.tif"), "wb") as handle:
            handle.write(b"junk")
        before = dict((n, digest(wat(n))) for n in os.listdir(wdir))
        code, out, err = run([wdir, "--json"])
        after = dict((n, digest(wat(n))) for n in os.listdir(wdir))
        check(before == after, "without --apply no file in the tree changes")
        code, out, err = run([wdir, "--apply", "--json"])
        doc = dict((os.path.basename(f["path"]), f) for f in
                   json.loads(out)["files"])
        check(doc["none.tif"]["verdict"] == OK
              and doc["none.tif"]["apply"] == "refreshed",
              "a raster with no statistics gets a sidecar and then passes")
        check(os.path.exists(wat("none.tif.aux.xml"))
              and not os.path.exists(wat("none.tif.aux.xml") + BAK_SUFFIX),
              "a new sidecar has nothing to back up")
        check(doc["stale.tif"]["verdict"] == OK, "a stale raster passes after")
        check(digest(old_side + BAK_SUFFIX) == old_bytes,
              "the old sidecar was copied byte for byte before the rewrite")
        new_text = digest(old_side).decode("utf-8")
        check("<Min>0</Min>" in new_text and 'key="K"' in new_text
              and "keep me" in new_text,
              "histograms, other domains and other keys are kept")
        check("STATISTICS_APPROXIMATE" not in new_text,
              "the approximate flag is dropped when exact numbers are written")
        check(new_text.count('key="STATISTICS_MAXIMUM"') == 1,
              "an existing key is updated, not repeated")
        check(digest(wat("stale.tif")) == raster_before,
              "--apply never touches the raster itself  <-- pinned defect")
        check([n for n in os.listdir(wdir) if n.endswith(".tmp")] == [],
              "no temporary file is left behind")
        check(doc["inside.tif"]["apply"].startswith("not refreshed: the NoData")
              and "NODATA_INSIDE_RANGE" in codes(doc["inside.tif"]),
              "a NoData value inside the range is not papered over"
              "  <-- pinned defect")
        check(not os.path.exists(wat("inside.tif.aux.xml")),
              "and no sidecar is written for it")
        check(doc["alln.tif"]["apply"] == "not refreshed: nothing to refresh",
              "an all NoData raster has no numbers to write")
        check("needs a decision" in doc["u8bad.tif"]["apply"],
              "a NoData value no cell can equal is not papered over")
        check(doc["clean.tif"]["apply"] == "not refreshed: nothing to refresh",
              "a clean raster is left alone")
        check(doc["broken.tif"]["apply"].endswith("could not be read"),
              "an unreadable raster is left alone")
        check(doc["nobandsidecar.tif"]["verdict"] == OK
              and 'domain="Q"' in digest(
                  wat("nobandsidecar.tif.aux.xml")).decode("utf-8"),
              "a sidecar with only a domain block gets a plain block beside it")
        check(doc["nobandfile.tif"]["verdict"] == OK
              and 'band="1"' in digest(
                  wat("nobandfile.tif.aux.xml")).decode("utf-8"),
              "a sidecar with no band gets band 1")
        side = wat("stale.tif.aux.xml")
        with open(side, "w") as handle:
            handle.write('<PAMDataset><PAMRasterBand band="1"><Metadata>'
                         '<MDI key="STATISTICS_MINIMUM">-5</MDI></Metadata>'
                         '</PAMRasterBand></PAMDataset>')
        stale_again = digest(side)
        code, out, err = run([wat("stale.tif"), "--apply"])
        check("a backup already exists" in out and code == 1
              and digest(side) == stale_again,
              "an existing backup is never overwritten")
        os.remove(side + BAK_SUFFIX)
        real_replace = os.replace

        def refuse(src, dst):
            raise OSError("disk full")
        os.replace = refuse
        try:
            code, out, err = run([wat("stale.tif"), "--apply"])
        finally:
            os.replace = real_replace
        check("could not write the sidecar: disk full" in out
              and digest(side) == stale_again,
              "a failed write leaves the old sidecar whole")
        check([n for n in os.listdir(wdir) if n.endswith(".tmp")] == [],
              "and removes its temporary file")

        # ---- text of the reports
        check(describe({"verdict": UNREADABLE, "path": "p", "reason": "r"})
              == ["UNREADABLE p  r"], "an unreadable report is one line")
        rep = {"verdict": WARN, "path": "p", "cells": 4, "valid_cells": 3,
               "nodata": None, "nodata_cells": 0,
               "findings": [{"code": "NO_STATS", "message": "m"}],
               "apply": "refreshed"}
        text = describe(rep)
        check("nodata=none" in text[0]
              and text[1].strip().startswith("NO_STATS")
              and text[2] == "    apply: refreshed",
              "a finding and an apply note are printed under the raster")

    # ---- the source itself
    with open(os.path.abspath(__file__), "rb") as handle:
        source = handle.read()
    check(all(b < 128 for b in bytearray(source)), "the source is pure ASCII")
    text = source.decode("ascii")
    banned = ("soc" + "ket", "url" + "lib", "http" + ".", "envi" + "ron",
              "subpro" + "cess")
    check(not any(word in text for word in banned),
          "the tool makes no network call, reads no variables, starts no program")
    check(not any(isinstance(node, ast.JoinedStr)
                  for node in ast.walk(ast.parse(text))),
          "no f-strings, for Python 3.9")

    # ---- the harness reports its own failures
    def boom():
        raise UnreadableRaster("boom")
    probe = _Tally()
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        probe.check(True, "a")
        probe.check(False, "b")
        probe.raises(boom, "c")
        probe.raises(lambda: 1 / 0, "d")
        probe.raises(lambda: None, "e")
        probe_code = probe.finish()
    check(probe_code == 1 and probe.passed == 2 and len(probe.failed) == 3,
          "a failed assertion is counted and makes the run exit 1")
    check("5 assertions, 3 failed" in sink.getvalue()
          and "FAILED: d (wrong exception" in sink.getvalue()
          and "FAIL  e (no error raised)" in sink.getvalue(),
          "the harness names each failed assertion and why")
    return tally.finish()


if __name__ == "__main__":  # pragma: no branch
    sys.exit(main())

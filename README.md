# nodatalie

Gate on stored versus recomputed raster statistics, and name the NoData value that explains the gap.

A 100 by 100 elevation raster has 500 cells that hold -9999, the old sentinel for "no measurement".
Someone runs the statistics once, with those cells counted: the minimum is -9999 and the mean is
-95.85.
Someone else then sets -9999 as the NoData value. That edit changes how every tool reads the cells.
It does not touch the stored statistics. GDAL, asked for the statistics of that file afterwards,
returns minimum -9999 and mean -95.85 again. The file is correct, the numbers are stale, and
nothing says so.

The opposite order is just as quiet. The statistics are stored with the sentinel masked, minimum
-50. Later the NoData value is cleared on a re-export. The 600 cells holding -9999 are now data.
GDAL still answers minimum -50, because that is what is stored. The real minimum is -9999 and the
real mean is 343 instead of 424. A colour ramp, a legend or an import check that reads the stored
numbers never looks at the cells, so it never sees the difference.

This tool reads the cells. It recomputes the minimum, maximum, mean and standard deviation with
the NoData cells masked, and compares them with the stored numbers. A gap is a finding and the
exit code is 1. It also says which gap it is: stored numbers that counted the NoData cells,
stored numbers that are simply old, or a NoData value that sits inside the real data range and
so masks real cells.

```
$ python nodatalie.py demo
FAIL       demo/dem_late.tif  cells=12 valid=10 nodata=-9999 nodata_cells=2
    STATS_COUNT_NODATA     stored statistics include the 2 cells equal to the declared NoData value -9999: either the NoData value changed after the statistics were stored, or it masks real cells
FAIL       demo/slope.asc  cells=12 valid=10 nodata=-9999 nodata_cells=2
    STALE_STATS            stored statistics differ from the cells: max stored 99, recomputed 18
WARN       demo/bathy.tif  cells=8 valid=5 nodata=0 nodata_cells=3
    NO_STATS               no stored statistics: a legend or stretch has nothing to read
    NODATA_INSIDE_RANGE    NoData value 0 lies inside the data range -5..40 and masks 3 cells (37.50 percent): confirm they are not real values
OK         demo/dem_ok.tif  cells=12 valid=10 nodata=-9999 nodata_cells=2
NODATALIE fail=2 warn=1 ok=1 unreadable=0 scanned=4
```

## What already exists

GDAL is the right tool for most raster questions, and this tool does not replace it. It reads
every format, compressed, tiled, multi-band and cloud-optimised files included, and it is much
faster. `gdalinfo -mm` forces a scan for the minimum and maximum. `gdalinfo -stats` prints the
statistics. Its documentation says it will "Force computation if no statistics are stored in an
image", and the source of `GDALRasterBand::GetStatistics` shows what that means: when
`STATISTICS_MINIMUM`, `STATISTICS_MAXIMUM`, `STATISTICS_MEAN` and `STATISTICS_STDDEV` are all
present, it returns them without reading a cell. The one exception is a stored set flagged
`STATISTICS_APPROXIMATE` when exact statistics were asked for. So the command people run to "check
the statistics" returns the stale numbers when they exist. The `rasterio` call `statistics()`
behaved the same way in the test below: it returned the stored minimum -50 for a raster whose
cells hold -9999.

Esri's Calculate Statistics tool recomputes, and its "Skip Existing" setting decides whether to
keep stored statistics or regenerate them. Its "Ignore Values" setting excludes values such as
NoData from the calculation.

What none of these does is compare. GDAL and Esri produce or print statistics. They do not read
the stored ones, recompute, diff, and exit non-zero. That comparison, run on a folder, is the one
thing this tool adds. A search of GitHub repositories for `nodatalie` returned no result on
2026-10-22, and searches for `raster statistics nodata validate`, `stale raster statistics
aux.xml` and `gdal check statistics nodata` returned no result either. That says nothing about
tools that are named differently.

## Requirements

Python 3.9 or newer. Standard library only: `struct`, `array`, `zlib`, `xml.etree.ElementTree`,
`argparse`, `json`, `math`. No GDAL, no `arcpy`, no `numpy`, no network, no install step.

```
git clone https://github.com/uhsear/nodatalie.git
```

## Quick start

```
python nodatalie.py --self-test
python nodatalie.py D:\rasters
```

To build the four files used above, save this as `make_demo.py` next to `nodatalie.py` and run
it. The helpers it calls are in `nodatalie.py` because the self-test uses them too.

```python
import os
import nodatalie as n

os.makedirs("demo", exist_ok=True)
dem = [[10, 12, 14, 16], [11, -9999, 15, 17], [12, 14, -9999, 18]]
flat = [v for row in dem for v in row]

n.make_tiff("demo/dem_ok.tif", dem, nodata="-9999")
n.make_sidecar("demo/dem_ok.tif", n.good_items(flat, -9999))

n.make_tiff("demo/dem_late.tif", dem, nodata="-9999")
n.make_sidecar("demo/dem_late.tif", n.good_items(flat))    # stored before NoData was set

n.make_tiff("demo/bathy.tif", [[-5, 0, 0, 9], [-3, 4, 0, 40]], nodata="0")   # no statistics

n.make_asc("demo/slope.asc", dem)
n.make_sidecar("demo/slope.asc", dict(n.good_items(flat, -9999), STATISTICS_MAXIMUM="99"))
```

## Usage

Run a check before any write. Nothing is written without `--apply`.

```
python nodatalie.py dem.tif
python nodatalie.py D:\rasters --strict
python nodatalie.py D:\rasters --json > report.json
python nodatalie.py D:\rasters --apply
```

| Flag | Default | What it does |
|---|---|---|
| `PATH ...` | none | `.tif`, `.tiff` or `.asc` files, or folders searched recursively. Required. |
| `--tol` | `1e-6` | Relative tolerance for each stored number. |
| `--require-stats` | off | A raster with no stored statistics fails instead of warning. |
| `--strict` | off | A warning exits 1 as well. |
| `--json` | off | One JSON document, worst raster first. |
| `--apply` | off | Write the recomputed statistics into the PAM sidecar. See below. |
| `--self-test` | off | Run the offline assertions and exit. |

Exit codes: 0 clean, 1 a finding (a warning too under `--strict`), 2 a raster could not be read,
a path does not exist or the folder holds no raster, 64 usage error. A finding outranks an
unreadable file. The parser sets `allow_abbrev=False`, so `--ap` is refused with exit 64 and a
short prefix can never turn on a write.

### Findings

| Code | Level | Meaning |
|---|---|---|
| `STALE_STATS` | FAIL | A stored number differs from the cells. The message names each one. |
| `STATS_COUNT_NODATA` | FAIL | The stored numbers match the cells only when the NoData cells are counted. The NoData value was set after the statistics, or it masks real cells. |
| `NO_STATS` | WARN | No stored statistics. FAIL under `--require-stats`. |
| `NODATA_INSIDE_RANGE` | WARN | The NoData value lies strictly between the smallest and largest real cell, and cells hold it. |
| `NODATA_UNREPRESENTABLE` | WARN | An integer raster declares a NoData value no cell can equal, such as -9999 on 8-bit data. |
| `ALL_NODATA` | FAIL | Every cell is NoData or NaN, so no statistic exists. |

A stale finding is only a WARN when the stored numbers carry `STATISTICS_APPROXIMATE`, because
approximate statistics are allowed to differ.

### What `--apply` writes

`--apply` writes one file per raster: `name.tif.aux.xml`, the GDAL PAM sidecar. It never opens a
raster for writing. It touches only rasters that failed on `STALE_STATS`, `STATS_COUNT_NODATA` or
`NO_STATS`, and refuses one that also has `NODATA_INSIDE_RANGE` or `NODATA_UNREPRESENTABLE`,
because fresh numbers would hide the NoData question. It refuses an `ALL_NODATA` raster too,
since there is nothing to write.

If a sidecar exists, `--apply` copies it to `name.tif.aux.xml.nodatalie.bak`, edits it in place
and keeps every other element (histograms, other metadata domains, other keys). It never
overwrites an existing backup: it stops and says so. The new file is written to a temporary name
and renamed, so a crash leaves one whole file. It then checks the raster again and prints the new
verdict. Format of the numbers: GDAL's own, 14 significant digits, and a valid percent rounded as
GDAL rounds it.

```
$ python nodatalie.py demo --apply
WARN       demo/bathy.tif  cells=8 valid=5 nodata=0 nodata_cells=3
    NO_STATS               no stored statistics: a legend or stretch has nothing to read
    NODATA_INSIDE_RANGE    NoData value 0 lies inside the data range -5..40 and masks 3 cells (37.50 percent): confirm they are not real values
    apply: not refreshed: the NoData value needs a decision first; refreshing would hide it
OK         demo/dem_late.tif  cells=12 valid=10 nodata=-9999 nodata_cells=2
    apply: refreshed
OK         demo/dem_ok.tif  cells=12 valid=10 nodata=-9999 nodata_cells=2
    apply: not refreshed: nothing to refresh
OK         demo/slope.asc  cells=12 valid=10 nodata=-9999 nodata_cells=2
    apply: refreshed
NODATALIE fail=0 warn=1 ok=3 unreadable=0 scanned=4
```

## How it decides

1. Read the cells in blocks and drop NaN. GDAL excludes NaN from statistics too, and the NaN
   cells never reach the count (observed, see Verification).
2. Mask the cells equal to the declared NoData value. A float32 raster compares in float32, so
   the text `-3.40282346638529e+38` still matches the cell `-3.4028234663852886e+38`.
3. Keep a running count, mean and sum of squared deviations per block and merge them. One pass,
   and no list of all cells in memory.
4. Take the standard deviation as the square root of the sum of squares divided by the count.
   GDAL divides by the count and not by count minus one. A recompute that divides by count minus
   one flags every healthy raster, so the self-test pins it.
5. Read the stored numbers: the in-file `GDAL_METADATA` tag for a TIFF, then the `.aux.xml`
   sidecar on top of it key by key. GDAL does the same, and the self-test pins it.
6. Compare the minimum, maximum, mean, standard deviation and valid percent that are stored. A
   missing key is not compared. Numbers match inside `--tol`, relative, with a floor of 1.0 on
   the scale so a mean near zero does not fail on noise.
7. On a gap, compare the stored numbers with the statistics of all cells, NoData counted. A match
   is `STATS_COUNT_NODATA`. No match is `STALE_STATS`.

## Verification

The self-test needs no data. It builds every raster itself, in pure Python, from a TIFF writer
and an LZW encoder that live in the same file.

```
$ python nodatalie.py --self-test
nodatalie self-test: no GDAL, no arcpy, no network, no data
--------------------------------------------------------------------
PASS  no values give the empty accumulator
PASS  count, min and max
PASS  mean and spread
PASS  an empty accumulator has no summary
PASS  the deviation divides by the count, as GDAL does  <-- pinned defect
PASS  the sample deviation is not what GDAL stores  <-- pinned defect
...
PASS  never rounds a partly masked raster up to 100  <-- pinned defect
...
PASS  a sample deviation is not accepted as the stored one  <-- pinned defect
...
PASS  a token split across read blocks is rejoined  <-- pinned defect
PASS  a last cell with no newline is read  <-- pinned defect
...
PASS  NoData added after the statistics were stored is named  <-- pinned defect
...
PASS  a read block smaller than one cell still advances  <-- pinned defect
PASS  LZW across the 9 to 12 bit code widths  <-- pinned defect
PASS  LZW on noise that fills the table and forces a reset  <-- pinned defect
...
PASS  a 15 digit float32 NoData text still matches the cells  <-- pinned defect
...
PASS  a sidecar overrides the tag, as GDAL does  <-- pinned defect
...
PASS  predictor 2 with no compression is refused  <-- pinned defect
...
PASS  a unique prefix --ap of --apply is REFUSED by the parser  <-- pinned defect
PASS  a longer prefix --app of --apply is REFUSED  <-- pinned defect
PASS  a unique prefix of --self-test is REFUSED  <-- pinned defect
PASS  a unique prefix of --json is REFUSED  <-- pinned defect
...
PASS  --apply never touches the raster itself  <-- pinned defect
...
PASS  a NoData value inside the range is not papered over  <-- pinned defect
...
--------------------------------------------------------------------
248 assertions, 0 failed
```

That run is 248 assertions on Windows under `C:\Python313\python.exe` 3.13.2. The same count
passes on Windows under 3.9.25, 3.10.20 and 3.12, and on Ubuntu under `python3` 3.12.3. The lines
cut with `...` are more of the same block above them.

Branch coverage of `nodatalie.py` is 100 percent: 1,326 statements and 354 branches, measured
with `coverage run --branch nodatalie.py --self-test`. One line carries `# pragma: no branch`,
the `if __name__ == "__main__":` guard. Coverage runs the file as a script, so the guard's false
arc can never be taken. `pyflakes` is clean and `python -W error` compiles the file clean.

Beyond the self-test, these were checked against GDAL itself (rasterio 1.5.1 on Windows), on
synthetic data:

| Claim | How it was exercised | Result |
|---|---|---|
| The TIFF reader and the fixture writer agree with GDAL cell for cell | 144 files: 8 sample types, both byte orders, no compression, LZW and deflate, predictor 2 on integers, strips and 16-pixel tiles. Each was read by GDAL and by this tool. | 0 differences |
| The recomputed statistics equal GDAL's | 7 files written by GDAL: LZW with predictor 2 and tiles, deflate, striped, float32 with NaN and a float32-maximum NoData, NaN as NoData, 8-bit with NoData 0 | min, max, mean and standard deviation equal to 1e-9 relative |
| GDAL returns stored statistics without scanning | A GDAL file got its statistics, then its NoData value was cleared. `statistics()` kept returning the old minimum. | min -50 returned, cells hold -9999 |
| The `STATS_COUNT_NODATA` case is real | A GDAL file got its statistics with no NoData value, then NoData -9999 was set. | GDAL still returned min -9999, mean -95.8476. This tool: `STATS_COUNT_NODATA`. |
| NaN is excluded from the statistics | float32 file of 6 cells: one NaN, one NoData | GDAL stored minimum 1 and valid percent 66.67, so 4 valid cells |
| A sidecar overrides the in-file tag, key by key | GDAL file with all four statistics in the tag and only `STATISTICS_MAXIMUM` in the sidecar | `tags()` showed the sidecar maximum and the three tag values |
| An ASCII grid reads the same as GDAL | 31 x 17 grid written by the fixture helper, read by both | the means differ by 6e-10 relative, the minimum is equal |
| Hostile sidecar XML is refused | A sidecar with a nested entity expansion | `UNREADABLE`, in 0.1 s, on Windows (expat 2.6.4) and Ubuntu (expat 2.6.1) |
| Speed | 2000 x 2000 int16, one run on the Windows workstation | 4.2 million cells a second uncompressed, 3.5 with deflate tiles, 0.7 with LZW and predictor, 2.4 for ASCII |

Every pinned defect marked in the self-test was reintroduced one at a time to confirm it fails
the run: the sample deviation, the valid-percent clamp, the token carry across read blocks, the
zero-size read block, the float32 cast, the tag and sidecar order, the abbreviated write flag,
the NoData-in-range block, the count-NoData branch, the LZW code width, the NaN count, the backup
guard, the two exit-code rules and the signed predictor. All 15 made the run fail, crash or hang.

## What it will not do

- Read most real-world rasters. It reads classic TIFF and Esri ASCII grid only. BigTIFF, JPEG,
  ZSTD and PackBits compression, the floating-point predictor, more than one band, Esri file
  geodatabase rasters, IMG, Esri GRID and cloud formats come back `UNREADABLE` with a reason.
  `gdalinfo -mm` is the answer for those.
- Match GDAL on speed. It is pure Python. A raster of 100 million cells takes about 25 seconds
  uncompressed and about 2.5 minutes with LZW and predictor 2.
- Know what the right NoData value is. `STATS_COUNT_NODATA` says the NoData value and the stored
  numbers disagree. A person decides which one is right.
- Read statistics stored anywhere but `GDAL_METADATA` and the `.aux.xml` sidecar. Statistics held
  by an Esri geodatabase or an `.ovr` file are invisible to it. It was not run against ArcGIS Pro
  output, because no Esri software was available.
- Read a NoData value from the sidecar. GDAL can keep one there for formats that cannot store it.
  This tool takes the value from the TIFF tag or the ASCII header only.
- Apply Esri's default. Esri documents `-9999` as the default `NODATA_value` of an ASCII grid
  when the header omits it. This tool reads a header with no `NODATA_value` as declaring none.
- Match GDAL to the last digit on ASCII grids. GDAL autodetects `float32` for many grids, and
  this tool reads them as `float64`. The difference is about 1e-8 relative, far inside the
  default tolerance.
- Defend against every hostile file. The sidecar is parsed with the standard library's
  `ElementTree`, which does not fetch external entities and stops runaway entity expansion on the
  expat versions listed above. An older expat is untested.
- Print strict JSON for every value. A stored value of `nan` in a sidecar comes out of `--json` as
  a bare `NaN`, which strict JSON parsers refuse. Python's `json` module reads it.

## Sources

- GDAL, `gdalinfo`, options `-stats`, `-mm`: https://gdal.org/en/stable/programs/gdalinfo.html
- GDAL raster data model, `STATISTICS_MINIMUM`, `STATISTICS_MAXIMUM`, `STATISTICS_MEAN`,
  `STATISTICS_STDDEV`, `STATISTICS_APPROXIMATE`, `STATISTICS_VALID_PERCENT`:
  https://gdal.org/en/stable/user/raster_data_model.html
- GDAL GeoTIFF driver, the `TIFFTAG_GDAL_NODATA` (42113) and `TIFFTAG_GDAL_METADATA` (42112) tags
  and the `.aux.xml` sidecar: https://gdal.org/en/stable/drivers/raster/gtiff.html
- GDAL `GDALRasterBand::GetStatistics`, `ComputeStatistics`, `SetStatistics` and
  `SetValidPercent`, read in `gcore/gdalrasterband.cpp` on the master branch on 2026-10-22. They
  show the stored values returned without a scan, the population standard deviation, the `%.14g`
  format and the valid-percent rounding: https://github.com/OSGeo/gdal/blob/master/gcore/gdalrasterband.cpp
- GDAL, `GDALRasterBand` class reference:
  https://gdal.org/en/stable/api/gdalrasterband_cpp.html
- GDAL ASCII grid driver: https://gdal.org/en/stable/drivers/raster/aaigrid.html
- Esri, ASCII raster format, the header keys and the optional `NODATA_value` with default -9999:
  https://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/esri-ascii-raster-format.htm
- Esri, Calculate Statistics, "Skip Existing" and "Ignore Values":
  https://doc.esri.com/en/arcgis-pro/latest/tool-reference/data-management/calculate-statistics.html

The `gdalinfo` command line was not installed on the test machines. Its behaviour here comes from
its documentation and from the GDAL source. The matching library call was run through `rasterio`.

## Contributing

Open an issue or pull request on GitHub.

## Author

Built by [Asir Khan](https://www.linkedin.com/in/asir-khan-310317264/).

## License

MIT.

## Related

Other single-file tools in this portfolio:

- [safe-republish](https://github.com/uhsear/safe-republish) - the same refusal before a destructive step, for a feature class replace
- [arcpy-nullscan](https://github.com/uhsear/arcpy-nullscan) - the placeholder-value problem in attribute tables, where -9999 hides too
- [gdbxray](https://github.com/uhsear/gdbxray) - schema a stock reader will not show you

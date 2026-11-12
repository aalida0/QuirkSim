# Recover the completed late32 transport without rerunning it

The reported run finished every transport case, then its audit process was
killed before local_audit.json was written. The signal source is not proven
without AF kernel/cgroup logs; excessive memory use is a plausible cause.
Source inspection found two avoidable full-history memory paths:

* Directory audits copied the complete ROOT file into BytesIO, then loaded
  every step column. ROOT handles/caches were not explicitly closed per case.
* Post-neutral packaging loaded full step histories again before finding a
  small return window.

The directory audit now opens ROOT files directly, disables object/array
caches, closes each file, and reduces endpoint checks in 8 MB chunks. Time
and proper-clock ordering are checked across chunk boundaries. Existing
macro/plan hashes, local decay conservation, shells, directions, handoffs,
completion and endpoint-count checks remain active.

The post-neutral packager recreates a missing audit from the existing files,
writing it atomically. A failed audit still prevents packaging. Its return
search streams neutral coordinates into temporary files, then merges the
two endpoint histories in bounded row blocks through disk-backed mappings.
This uses the original piecewise-linear search, not downsampling. A second
ROOT pass loads only the selected +/-0.5 ps window. Temporary coordinate
files are removed after each case. The final compact archive retains the
population/censoring inventory and excludes ROOT files.

Long-duration automatic packaging and --pack-only route to this path.
Neither recovery command below invokes execute(), compilation or transport:

```bash
cd ~/QuirkSim
git pull --ff-only origin research/quirkpheno-weak-recoil-pilot
python scripts/package_triplet_postneutral.py \
  /data/aalida/QuirkSim/triplet_late32_20261112T001756Z
```

Progress prints Auditing ... and then Packing ... . Return the printed
*_postneutral.tar.gz file. The existing ROOT files remain in place.

Validation includes one-row audit chunks versus the original whole-test
reference, asynchronous endpoint grids across seven-row search blocks,
exact selected-window comparisons, ROOT-to-compact round trips, retained
censored counts and recreation of the missing-audit control path. These
are local synthetic/regression checks; the large AF archive has not been
read here and successful AF packaging is not claimed in advance.

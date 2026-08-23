Server build of `probabilisticNeighbour`, for the bare-metal box (no scheduler).

The physics is identical to `src/probabilisticNeighbour/probabilisticNeighbour.cpp`. Only
the plumbing differs: the output path is the server scratch directory, the snapshot count
is a command-line argument (snapshots dominate disk at large `L`), and the largest-spot
series is sampled every 100 sweeps instead of every 10.

## Running it

```bash
scp -r probabilisticNeighbourServer/ server:~/            # or clone the repo
cd probabilisticNeighbourServer
DRY_RUN=1 ./run_server.sh        # prints the plan, runs nothing
nohup ./run_server.sh &> run.out &
```

The script compiles the binary itself if it is missing or out of date, so nothing needs
building by hand. Output goes to `/nbi/home/rpw391/cell-disk/solarFlares/outputs`; override with
`OUTDIR=... ./run_server.sh` for a test run somewhere else.

**It resumes.** A run counts as finished when its `spotSize_*.tsv` exists and is
non-empty, and that file is written only after the sweep loop completes — so a killed or
crashed run is never mistaken for a good one. If the machine goes down, just start the
script again and it will pick up what is missing.

## The grid

| L | rho | p | seeds | runs | cost |
|---|---|---|---|---|---|
| 128 | 0.2, 0.4, 0.6, 0.8 | 0:0.05:1 | 64 | 5376 | 0.97M core-s |
| 256 | 0.2, 0.6 | 0:0.05:1 | 32 | 1344 | 0.86M core-s |
| 512 | 0.2 | 0:0.05:1 | 24 | 504 | 1.42M core-s |
| 1024 | 0.2 | 11 values | 16 | 176 | 1.81M core-s |

~7400 runs, ~5.1M core-seconds, about **11 h on 124 cores**, longest single job under 3 h.
`SEED_MULT=2` roughly doubles all of it.

## Disk and transfer

Every output is gzipped as soon as its run finishes, so the uncompressed set never exists
on disk. The sign lattice is mostly zeros and large domains, so it packs well: a measured
22 MB snapshot file goes to 1.7 MB (12x). Histograms compress ~3x, the largest-spot series
~2x. Snapshot counts are also scaled down per L (4 at L=128 rising to 8 at L=512), since
64 seeds x 4 snapshots is already 256 independent configurations per (rho, p).

Together that puts the whole sweep at roughly **1 GB**, against ~11 GB uncompressed.

Measured cost per 1e5 sweeps (laptop, p=0.9): 90 s at L=128, 320 s at L=256, 1410 s at
L=512, 5150 s at L=1024, and **54,500 s at L=2048**. The jump at the top is not L^2 — the
lattice stops fitting in cache and the random access pattern starts missing — which is why
L=2048 is excluded by default. One such run would take ~30 h on the server on its own.
`BIG=1` enables it if you ever want it.

Cost also depends on `p`: the neighbour rule scans the occupied list to find a like-sign
spot with an empty neighbour, so p=1 is dearer than p=0. The figures above use p=0.9, so
the total is a conservative estimate.

## Getting the data back

Filenames match the laptop convention, and `plots.py` now reads `.tsv.gz` as happily as
`.tsv`, so nothing needs unpacking:

```bash
rsync -av server:/nbi/home/rpw391/cell-disk/solarFlares/outputs/ src/probabilisticNeighbour/outputs/
cd src/probabilisticNeighbour && python plots.py
```

Leave the files gzipped -- `numpy.loadtxt` decompresses transparently and the snapshot
reader uses `gzip.open`. If you ever do want them expanded, `gunzip *.gz` in `outputs/`.

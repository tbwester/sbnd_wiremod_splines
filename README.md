# SBND WireMod Splines

Calculate N-dimensional splines for SBND WireMod

## Building

Build within an SL7 container. Run

```
source build_env.sh
mkdir build && cd build
cmake /path/to/sbnd_wiremod_splines
cmake --build .
```

## Grid submission

Do not run in an SL7 container. First, set up the fife_utils package Spack:

```
source "/cvmfs/fermilab.opensciencegrid.org/packages/common/setup-env.sh"
spack load /bkga3gp # fife-utils
```

Modify the output paths and data set in `grid/job.cfg`, then run

```
fife_launch -c /path/to/job.cfg --stage=hist
```

## Output files

Depending on the number of jobs indicated in `job.cfg`, output histograms will be split into multiple files. With ROOT loaded, the histograms can be combined into a single file via

```
hadd -f output.root /path/to/job/output/*.root
```

## Analysis

See the `macros/` directory for creating projections and fitting the N-dimensional histograms.

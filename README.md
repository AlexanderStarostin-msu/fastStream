# fastStream
Trajectory analysis in the flow field

## Building `fslFoam`

In addition to OpenFOAM 10, the application requires Boost.Geometry headers.
On Ubuntu, install them with:

```sh
sudo apt-get install libboost-dev
```

Boost.Geometry is header-only; no Boost library needs to be linked separately.
Then source the OpenFOAM environment and build:

```sh
source /path/to/OpenFOAM-10/etc/bashrc
cd src/fslFoam
wmake
```

The shared `minObstacleClearance` setting in `fsStreamLineDict` defaults to
zero, which is valid and disables obstacle expansion. The retired
`deviceRadius` setting is not read by the current application.

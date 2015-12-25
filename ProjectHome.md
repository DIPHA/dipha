# DIPHA (A Distributed Persistent Homology Algorithm), v2.1.0 #
Copyright 2014 IST Austria
## Project Founder: ##

Jan Reininghaus (Email: jan.reininghaus@gmail.com)


## Contributors: ##

Ulrich Bauer, Michael Kerber

## Downloads: ##

  * [DIPHA, v2.1.0](https://docs.google.com/uc?id=0B7Yz6TPEpiGEYTZiMHkxZVpLUU0&export=download)

## Description: ##

This C++ software package computes persistent homology in a distributed setting following the algorithm proposed in `[`[1](http://dx.doi.org/10.1137/1.9781611973198.4)`]`. For an introduction to persistent homology, see the textbook `[`[2](http://www.ams.org/bookstore-getitem/item=mbk-69)`]`.

There are three types of input that are currently supported by `DIPHA`:

  1. d-dimensional gray-scale image data. The data is internally interpreted as a weighted cubical cell complex and its lower-star filtration is used for the subsequent persistence computation. This approach is described in detail in `[`[3](http://link.springer.com/chapter/10.1007%2F978-3-642-23175-9_7)`]`.
  1. distance matrix data. The data is internally interpreted as a Rips complex of as many points as there are columns in the given matrix. The distance between any two points is then defined by the input data.
  1. weighted regular cell complexes in (co-)boundary matrix form. This (fallback) input type allows the computation of persistent homology of e.g. alpha shapes, rips complexes, or witness complexes. The user is responsible for making sure that the weights of the cells induce a filtration of the complex.
The output produced by `DIPHA` consists of the persistence diagram. Each point in the diagram is defined by the dimension of the homological feature that it represents and the corresponding birth- and death value.

Input and output is realized using binary files whose format is specified below. `DIPHA` includes `MATLAB` functions to create the input files and visualize the output.

To achieve good performance `DIPHA` supports dualized computation as described in `[`[4](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.225.5421)`]`, makes use of the optimization introduced in `[`[5](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.224.6560)`]` and employs an efficient data structure developed in the [PHAT](http://phat.googlecode.com/) project.

## Setup: ##

Prerequisites:
  * A modern C++ compiler like [GCC](http://gcc.gnu.org/) 4.7, [Clang](http://clang.llvm.org/) 3.3, or [Visual Studio](http://www.microsoft.com/en-us/download/details.aspx?id=34673) 2012. `DIPHA` uses some C++11 features, so older compilers are not supported. If you are using Visual Studio 2013 you need to additionally install Visual Studio 2012.
  * An `MPI` implementation like [Open MPI](http://www.open-mpi.org/)  or [MPICH2](http://www.mpich.org/). Other MPI implementations may also work, but are untested.
  * the cross-platform, open-source build system [CMake](http://www.cmake.org/), version 2.8 or later

To compile `DIPHA`:
  1. use CMake to create a build environment. If you are using Visual Studio 2013 you need to run CMake twice.
  1. compile `DIPHA` using your favorite C++ compiler.

## Usage: ##

To run `DIPHA` using a single process:
```
dipha [options] input_filename output_filename
```
where the available `options` are:
  * `--benchmark`: prints some profiling information.
  * `--dual`: runs the dualized version of the algorithm, which dramatically improves the running time in some cases.
  * `--upper_dim D`: restricts the computation to the first `D` dimensions of the input complex. This option is mandatory when dealing with distance matrix data.
  * `--upper_value X`: computes the persistence diagram only up to value `X`. This option is crucial when computing persistence of large distance matrices.

To run `DIPHA` with N processes using `MPI`:
```
mpiexec -n N dipha [options] input_filename output_filename
```

## File Formats: ##

All file formats in `DIPHA` are in little-endian binary format. Integral values are stored using the binary representation of 64bit signed integers, while floating point values are stored in IEEE 754 double-precision binary floating-point format.

The first symbol in every `DIPHA` file is the magic number 8067171840. The second symbol is the integer encoding the actual file type. See `dipha/file_types.h` for a list of file types and their associated integral identifier. The rest of the symbols in the file depend on the file type:

  * d-dimensional gray-scale image data (`IMAGE_DATA`):
    1. number of data values `n`
    1. dimension `d`
    1. lattice resolution: `g`<sub>1</sub> ... `g`<sub>d</sub>
    1. floating point data values in x-fastest order: `v`<sub>1</sub> ... `v`<sub>n</sub>
> > Example: a gray-scale FullHD frame would have `n = 1920 * 1080`, `d = 2`, `g`<sub>1</sub> `= 1920`, `g`<sub>2</sub> `= 1080`, and `v`<sub>1</sub> ... `v`<sub>n</sub> would be the gray-scale values in scan-line order.
> > A simple `MATLAB` file writer is included in `DIPHA` (`matlab/save_image_data.m`).

  * weighted regular cell complexes in (co-)boundary matrix form (`WEIGHTED_BOUNDARY_MATRIX`):
    1. 0 if the file contains the boundaries of the cells, 1 if the file contains the coboundaries
    1. total number of cells `n`
    1. dimension of the cell complex `d`
    1. dimensions of the cells: `dim`<sub>1</sub> ... `dim`<sub>n</sub>
    1. floating point values of the cells: `value`<sub>1</sub> ... `value`<sub>n</sub>
    1. offsets of (co-)boundary entries of the cells: `offset`<sub>1</sub> ... `offset`<sub>n</sub>. Each offset represents the first index of the (co-)boundary in the flattened (co-)boundary matrix.
    1. number of entries in the (co-)boundary matrix m
    1. entries of the flattened (co-)boundary matrix m: `entry`<sub>1</sub> ... `entry`<sub>m</sub>
> > To convert a file from boundary format (required for primal computation) to coboundary format (required for dual computation), you may use the `dualize` utility.

  * distance matrix data (`DISTANCE_MATRIX`)
    1. number of input points `n`
    1. floating point values for the distances: `d`<sub>1,1</sub> ... `d`<sub>1,n</sub> `d`<sub>2,n</sub> ... `d`<sub>n,n</sub>
> > A simple `MATLAB` file writer is included in `DIPHA` (`matlab/save_distance_matrix.m`).

  * persistence diagram (`PERSISTENCE_DIAGRAM`):
    1. number of points in the diagram `p`
    1. dimension, birth, death triples: `dim`<sub>1</sub> `birth`<sub>1</sub> `death`<sub>1</sub> ... `dim`<sub>p</sub> `birth`<sub>p</sub> `death`<sub>p</sub>. Dimensions are stored as integers, while birth and death values are stored as floating point numbers.
> > If `dim`<sub>k</sub> is negative, then the corresponding triple represents an essential class of dimension -`dim`<sub>k</sub> - 1.
> > A simple `MATLAB` function to plot persistence diagrams is included in `DIPHA` (`matlab/plot_persistence_diagram.m`).

## Examples: ##
To save space, `DIPHA` does not come with any example data files. However, there are several `MATLAB` functions (`matlab/create_*.m`) that may be used to create synthetic examples for benchmark purposes. For example, to create a 3-dimensional image data set with 32<sup>3</sup> voxels whose values are given by random numbers you can execute the command
```
matlab -nojvm -nodisplay -nosplash -r "addpath('../matlab'); create_noise_image_data(3, 32); exit"
```
in the `dipha/examples` folder.

To visualize persistence diagrams with a large number of points it is recommended to use the provided `MATLAB` function that draws a density estimate of the diagram (`matlab/plot_persistence_diagram_density.m`).

## References: ##
`[`[1](http://dx.doi.org/10.1137/1.9781611973198.4)`]` U.Bauer, M.Kerber, J.Reininghaus: _Distributed computation of persistent homology_. Proceedings of the Sixteenth Workshop on Algorithm Engineering and Experiments (ALENEX), 2014.

`[`[2](http://www.ams.org/bookstore-getitem/item=mbk-69)`]` H.Edelsbrunner, J.Harer: _Computational Topology, An Introduction_. American Mathematical Society, 2010, ISBN 0-8218-4925-5.

`[`[3](http://link.springer.com/chapter/10.1007%2F978-3-642-23175-9_7)`]` H.Wagner, C.Chen, E.Vucini: _Efficient Computation of Persistent Homology for Cubical Data_. Topological Methods in Data Analysis and Visualization II, 2012.

`[`[4](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.225.5421)`]` V.de Silva, D.Morozov, M.Vejdemo-Johansson: _Dualities in persistent (co)homology_. Inverse Problems 27, 2011.

`[`[5](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.224.6560)`]` C.Chen, M.Kerber: _Persistent Homology Computation With a Twist_. 27th European Workshop on Computational Geometry, 2011.

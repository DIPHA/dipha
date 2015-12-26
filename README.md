# DIPHA (A Distributed Persistent Homology Algorithm)
Copyright 2014 IST Austria
## Project Founder: 

Jan Reininghaus (Email: jan.reininghaus@gmail.com)


## Contributors: 

Ulrich Bauer, Michael Kerber


## Description:

This C++ software package computes persistent homology following the algorithm proposed in [[A]](http://dx.doi.org/10.1137/1.9781611973198.4). Besides supporting parallel execution on a single machine, DIPHA may also be run on a cluster of several machines using MPI. For an introduction to persistent homology, see the textbook [[B]](http://www.ams.org/bookstore-getitem/item=mbk-69). 

There are three types of input that are currently supported by DIPHA:	

  1. d-dimensional gray-scale image data. The data is internally interpreted as a weighted cubical cell complex and its lower-star filtration is used for the subsequent persistence computation. This approach is described in detail in [[C]](http://link.springer.com/chapter/10.1007%2F978-3-642-23175-9_7).
  
  2. distance matrix data. The data is internally interpreted as a [Vietoris–Rips complex](https://en.wikipedia.org/wiki/Vietoris%E2%80%93Rips_complex) of as many points as there are columns in the given matrix. The distance between any two points is then defined by the input data. To allow for large data analysis, DIPHA supports sparse distance matrices which only encode the finite distances to save space.

  3. weighted regular cell complexes in (co-)boundary matrix form. This (fallback) input type allows the computation of persistent homology of e.g. Alpha complexes, Vietoris–Rips complex approximations, or Witness complexes. 
  
The output produced by DIPHA consists of the persistence diagram. Each point in the diagram is defined by the dimension of the homological feature that it represents and the corresponding birth- and death value. 

Input and output is realized using binary files whose format is specified below. DIPHA includes MATLAB functions to create the input files and visualize the output.

To achieve good performance DIPHA supports dualized computation as described in [[D]](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.225.5421), makes use of the optimization introduced in [[E]](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.224.6560) and employs an efficient data structure developed in the [PHAT project](http://phat.googlecode.com/) as described in [[F]](https://people.mpi-inf.mpg.de/~mkerber/bkrw-pphat.pdf).

## Setup:

Prerequisites:

  * A modern C++ compiler like [GCC 4.7](http://gcc.gnu.org/), [Clang 3.3](http://clang.llvm.org/), or [Visual Studio 2015](https://www.visualstudio.com/vs-2015-product-editions). DIPHA uses C++11, so older C++ compilers are not supported. In particular, Visual Studio 2013 is not supported.

  * An `MPI` implementation like [Open MPI](http://www.open-mpi.org/), [MPICH2](http://www.mpich.org/), or [MS-MPI](https://msdn.microsoft.com/en-us/library/vs/alm/bb524831(v=vs.85).aspx). Other MPI implementations may also work, but are untested.

  * the cross-platform, open-source build system [CMake](http://www.cmake.org/), version 2.8 or later

To compile DIPHA:

  1. use CMake to create a build environment.
  
  2. compile DIPHA using your favorite C++ compiler.

## Usage:

To run DIPHA using a single process:
```
dipha [options] input_filename output_filename
```
where the available `options` are:

  * `--benchmark`: prints some profiling information.

  * `--dual`: runs the dualized version of the algorithm, which dramatically improves the running time in some cases.
  
  * `--upper_dim D`: restricts the computation to the first `D` dimensions of the input complex. This option is mandatory when dealing with distance matrix data.

To run DIPHA with N processes using `MPI`:
```
mpiexec -n N dipha [options] input_filename output_filename
```

## File Formats:

All file formats in DIPHA are in little-endian binary format. Integral values are stored using the binary representation of 64bit signed integers, while floating point values are stored in IEEE 754 double-precision binary floating-point format.

The first symbol in every DIPHA file is the magic number 8067171840. The second symbol is the integer encoding the actual file type. See `dipha/file_types.h` for a list of file types and their associated integral identifier. The rest of the symbols in the file depend on the file type:

  * d-dimensional gray-scale image data (`IMAGE_DATA`):
    1. number of data values `n` 
    2. dimension `d` 
    3. lattice resolution: `g(1) ... g(d)` 
    4. floating point data values in x-fastest order: `v(1) ... v(n)` 
	
    Example: a gray-scale FullHD frame would have `n = 1920 * 1080`, `d = 2`, `g(1) = 1920`, `g(2) = 1080`, and `v(1) ... v(n)`, would be the gray-scale values in scan-line order.
    A simple MATLAB file writer is included in DIPHA (`matlab/save_image_data.m`).
	
  * weighted regular cell complexes in (co-)boundary matrix form (`WEIGHTED_BOUNDARY_MATRIX`):
    1. 0 if the file contains the boundaries of the cells, 1 if the file contains the coboundaries 
    2. total number of cells `n`
    3. dimension of the cell complex `d`
    4. dimensions of the cells: `dim(1) ... dim(n)`
    5. floating point values of the cells: `value(1) ... value(n)`
    6. offsets of (co-)boundary entries of the cells: `offset(1) ... offset(n)`. Each offset represents the first index of the (co-)boundary in the flattened (co-)boundary matrix.
    7. total number of non-zero entries `m` of the (co-)boundary matrix 
    8. entries of the flattened (co-)boundary matrix: `entry(1) ... entry(m)`
	
  To convert a file from boundary format (required for primal computation) to coboundary format (required for dual computation), you may use the included `dualize` utility.
    
  * distance matrix data (`DISTANCE_MATRIX`)
    1. number of input points `n` 
    2. floating point values for the distances in row order: `d(1,1) ... d(1,n) d(2,1) ... d(n,n)`
	
    A simple MATLAB file writer is included in DIPHA (`matlab/save_distance_matrix.m`).

  * sparse distance matrix data (`SPARSE_DISTANCE_MATRIX`)
    1. number of input points `n`
    2. number of connections of each point: `c(1) ... c(n)`
    3. indices of connected vertices: `i(1,1) ... i(1,c(1)) i(2,1) ... i(n,c(n))`
    4. distance values of connected vertices `d(1,1) ... d(1,c(1) d(2,1) ... d(n,c(n))`
    
 To create a sparse distance matrix DIPHA provides a `full_to_sparse_distance_matrix` utility.

  * persistence diagram (`PERSISTENCE_DIAGRAM`):
    1. number of points in the diagram `p`
    2. dimension, birth, death triples: `dim(1) birth(1) death(1) ... dim(p) birth(p) death(p)`. Dimensions are stored as integers, while birth and death values are stored as floating point numbers.
	
    If `dim(k)` is negative, then the corresponding triple represents an essential class of dimension `-dim(k) - 1`.
    A simple MATLAB function to plot persistence diagrams is included in DIPHA (`matlab/plot_persistence_diagram.m`).
    
## Examples:
To save space, DIPHA does not come with any example data files. However, there are several MATLAB functions (`matlab/create_*.m`) that may be used to create synthetic examples for benchmark purposes. For example, to create a 3-dimensional image data set with 32^3 voxels whose values are given by random numbers you can execute the command
```
matlab -nojvm -nodisplay -nosplash -r "addpath('../matlab'); create_noise_image_data(3, 32); exit"
```
in the `dipha/examples` folder.

To visualize persistence diagrams with a large number of points it is recommended to use the provided MATLAB function that draws a density estimate of the diagram (`matlab/plot_persistence_diagram_density.m`).
 
## References:
*A*) U.Bauer, M.Kerber, J.Reininghaus: [Distributed computation of persistent homology](http://dx.doi.org/10.1137/1.9781611973198.4). Proceedings of the Sixteenth Workshop on Algorithm Engineering and Experiments (ALENEX), 2014.

*B*) H.Edelsbrunner, J.Harer: [Computational Topology, An Introduction](http://www.ams.org/bookstore-getitem/item=mbk-69). American Mathematical Society, 2010, ISBN 0-8218-4925-5.

*C*) H.Wagner, C.Chen, E.Vucini: [Efficient Computation of Persistent Homology for Cubical Data](http://link.springer.com/chapter/10.1007%2F978-3-642-23175-9_7). Topological Methods in Data Analysis and Visualization II, 2012.

*D*) V.de Silva, D.Morozov, M.Vejdemo-Johansson: [Dualities in persistent (co)homology](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.225.5421). Inverse Problems 27, 2011.

*E*) C.Chen, M.Kerber: [Persistent Homology Computation With a Twist](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.224.6560). 27th European Workshop on Computational Geometry, 2011.

*F*) U.Bauer, M.Kerber, J.Reininghaus, H.Wagner: [PHAT – Persistent Homology Algorithms Toolbox](https://people.mpi-inf.mpg.de/~mkerber/bkrw-pphat.pdf). Mathematical Software – ICMS 2014, Lecture Notes in Computer Science Volume 8592, 2014, pp 137-143
  

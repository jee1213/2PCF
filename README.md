# kdtree 
 The set of tree codes here are implemented to compute the number of particles 
 in a given separation for all the particles in a 3-dimensional particle distribution data. 
 It follows the KD tree algorithm introduced in Moore, A., Connolly, A., Genovese, C. et al. (2000),
 "Fast Algorithms and Efficient Statistics: N-point Correlation Functions".

## Usage
 One is expected to have two sets of input particle distributions in 3d: (clustering) data and random. 
 The code can computes auto pair counts (data-data or random-ramndom) and cross pair counts (data-random)
 in a given (uniform) binning in 1d (radial) or 2d (radial, azimuthal).  
 The output can then be used to compute the two-point correlation function (2PCF) using Landy-Szalay (LS) estimator
 such that the clustering signal of the data can be estimated in comparison to the random particle distributions.
 An example IDL code that computes 2PCF using LS estimator is present (.pro).

### Requirements
 The geometry of the random and the clustering data is expected to be identical.
 tree_polar.c and tree_rr.c computes the number of pairs in two-dimensional bins of polar coordinates (r, theta)
 while tree_bin.c and tree_radial computes the number of pairs in one-dimensional bins in radial direction (r).
 
## Acknowledgement
 The initial code was developed by Nicolas Canac in python. Inh Jee has implemented it in C for the purpose of
 memory management for a large number of particles (on the orders of ~ millions), extensively debugged for a reliable
 performance, and parallelized it using OpenMP for a faster computation. 

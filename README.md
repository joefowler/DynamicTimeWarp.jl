# DynamicTimeWarp

[![Build Status](https://travis-ci.org/joefowler/DynamicTimeWarp.jl.svg?branch=master)](https://travis-ci.org/joefowler/DynamicTimeWarp.jl)

Dynamic Time Warping

A method for aligning sequences in a way that is insensitive to "warping" (stretching
or shifting) along the time axis.


To Do List:
* DTW Barycenter Averaging (see paper by Petitjean et al. 2011).
* Restricted DTW, where only a parallelogram is explored, not the full rectangular
  space. (This is a speed/memory optimization.)
* FastDTW (Salvador & Chan, Intelligent Data Analysis 2007), which works coarsening until
  DTW is easy, and then stepping back to the finest scales one at a time in a way that
  lets us search only "a radius" around the best DTW path from the next coarser scale.
* Possibly a "spectrum aligner" where the data being warped are histograms. (Does this
  affect anything? Maybe not needed.)
* Learn the terminology and actually get proper references in place!
* Some demonstration examples.
  
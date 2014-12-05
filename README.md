# DynamicTimeWarp

[![Build Status](https://travis-ci.org/joefowler/DynamicTimeWarp.jl.svg?branch=master)](https://travis-ci.org/joefowler/DynamicTimeWarp.jl)

Dynamic Time Warping

A method for aligning sequences in a way that is insensitive to "warping" (stretching
or shifting) along the time axis.


To Do List:
* DTW Barycenter Averaging (see paper by Petitjean et al. 2011).
* Restricted DTW, where only a parallelogram is explored, not the full rectangular
  space. (This is a speed/memory optimization.)
* Possibly a "spectrum aligner" where the data being warped are histograms. (Does this
  affect anything? Maybe not needed.)
* Learn the terminology and actually get proper references in place!
* Some demonstration examples.
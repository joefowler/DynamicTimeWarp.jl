# DynamicTimeWarp

[![Build Status](https://travis-ci.org/joefowler/DynamicTimeWarp.jl.svg?branch=master)](https://travis-ci.org/joefowler/DynamicTimeWarp.jl)

## Dynamic Time Warping

A method for aligning sequences in a way that is insensitive to
arbitrary "warping" (stretching or shifting) along the time (or "time") axis.
[Wikipedia article](http://en.wikipedia.org/wiki/Dynamic_time_warping)


## Features

* Standard DTW. This algorithm requires O(MN) time and memory, where M, N are the lengths
  of the two sequences to be aligned.
* FastDTW (Salvador & Chan, Intelligent Data Analysis 2007), which works by
  coarsening until DTW is easy, and then stepping back to the finest scales one at a
  time in a way that
  lets us search only "a radius" around the best DTW path from the next coarser scale.


## To Do List

- [*] Restricted DTW, where only a parallelogram is explored, not the full rectangular
  space. (This is a speed/memory optimization.)
- [*] FastDTW (Salvador & Chan, Intelligent Data Analysis 2007), which works by
  coarsening until DTW is easy, and then stepping back to the finest scales one at a
  time in a way that
  lets us search only "a radius" around the best DTW path from the next coarser scale.
- [*] Improve FastDTW code by creating a WindowedMatrix type, so that element access can
   be easier in the dtwwindowed code.
- [ ] Consider implementing SparseDTW [paper on arXiv](http://arxiv.org/abs/1201.2969)
- [ ] DTW Barycenter Averaging (see paper by Petitjean et al. 2011).
- [ ] Possibly a "spectrum aligner" where the data being warped are histograms. (Does this
  affect anything? Maybe not needed.)
- [ ] Learn the terminology and actually get proper references in place!
- [ ] Some demonstration examples.
  
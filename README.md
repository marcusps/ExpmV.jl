# ExpmV

[![Build Status](https://travis-ci.org/matteoacrossi/ExpmV.jl.svg?branch=master)](https://travis-ci.org/matteoacrossi/ExpmV.jl)

> This is a fork from https://github.com/marcusps/ExpmV.jl, implementing the evaluation for multiple values of the parameter ` t`.

This is a Julia translation of the MATLAB implementation of Al-Mohy and Higham's
function for computing `expm(t*A)*v` when `A` is sparse, without explicitly computing `expm(A)`.

If `t` is a `StepRangeLen` object (i. e. a `linspace`), use an optimized algorithm to output the result for all `t`.

The original code can be found at (https://github.com/higham/expmv), and the theory is explained in the following article:

*Computing the Action of the Matrix Exponential, with an Application to Exponential Integrators*, Awad H. Al-Mohy and Nicholas J. Higham, SIAM Journal on Scientific Computing 2011 33:2, 488-511. ([preprint](http://eprints.ma.man.ac.uk/1426/))


## Installation

Install into Julia using the package manager:

```julia
Pkg.add("ExpmV")
```

 ## Usage

```julia
expmv(t, A, v)
```

Eg. `t = 1.`, or `t = linspace(0, 1, 100)`.

## Benchmarks

Need to run updated benchmarks on Julia v1.0


## License

Released under the [BSD 2-clause license](https://tldrlegal.com/license/bsd-2-clause-license-(freebsd)) used in Al-Mohy and  Higham's original code, with the exception of function `norm1est`: the fast 1-norm estimation that is crucial for the speed of this algorithm is adapted from code available in Julia's `Base` module. This function (`norm1est`) is licensed under the MIT license.

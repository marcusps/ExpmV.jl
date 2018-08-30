# ExpmV

[![Build Status](https://travis-ci.org/matteoacrossi/ExpmV.jl.svg?branch=master)](https://travis-ci.org/matteoacrossi/ExpmV.jl)

> This is a fork from https://github.com/marcusps/ExpmV.jl, implementing the evaluation for multiple values of the parameter ` t`.

This is a Julia translation of the MATLAB implementation of Al-Mohy and Higham's
function for computing `expm(t*A)*v` when `A` is sparse, without explicitly computing `expm(A)`.

If `t` is a `StepRangeLen` object (i. e. a `linspace`), use an optimized algorithm to output the result for all `t`.

The original code can be found at (https://github.com/higham/expmv), and the theory is explained in the following article:

*Computing the Action of the Matrix Exponential, with an Application to Exponential Integrators*, Awad H. Al-Mohy and Nicholas J. Higham, SIAM Journal on Scientific Computing 2011 33:2, 488-511. ([preprint](http://eprints.ma.man.ac.uk/1426/))

The fast 1-norm estimation that is crucial for the speed of this algorithm is adapted from code available in Julia's `Base` module. This function (`norm1est`) is licensed under the MIT license.

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

This benchmark shows the performance of `ExpmV` compared to [Expokit.jl](https://github.com/acroy/Expokit.jl) and the builtin dense `expm` of Julia, for complex, non-Hermitian matrices. The benchmark is done using [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) on a Macbook Pro 2016 with 2,9 GHz Intel Core i5 and 16 GB RAM. The script is in `benchmark/compare.jl`.

### Matrix density 0.01
| Matrix rows                   |  `Expm` | `Expokit` | `Expmv`               |
|-----------------------|---------:|-----------:|-------:|
| 32      | 158.495 μs  |  30.100 μs  |  53.609 μs|
| 64      | 856.923 μs  |  52.036 μs  |  58.536 μs|
| 128     |   7.805 ms  | 537.083 μs  |  80.650 μs|
| 256     |  40.027 ms  |   2.993 ms  | 112.047 μs|
| 512     | 277.680 ms  |   3.195 ms  | 218.490 μs|
| 1024    |    1.902 s  |   4.267 ms  | 571.590 μs|

### Matrix density 0.001
| Matrix rows                   |  `Expm` | `Expokit` | `Expmv`               |
|-----------------------|---------:|-----------:|-------:|
| 32      |  31.147 μs  | 12.144 μs  | 55.103 μs |
| 64      | 471.424 μs  | 15.816 μs  | 53.599 μs |
| 128     |   7.368 ms  | 34.339 μs  | 60.320 μs |
| 256     |  27.817 ms  | 61.137 μs  | 76.773 μs |
| 512     | 325.282 ms  |182.016 μs  |142.402 μs |
| 1024    |    1.568 s  |  2.137 ms  |306.293 μs |

## License

Released under the [BSD 2-clause license](https://tldrlegal.com/license/bsd-2-clause-license-(freebsd)) used in Al-Mohy and  Higham's original code.

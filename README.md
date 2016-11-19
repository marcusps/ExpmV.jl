# ExpmV

[![Build Status](https://travis-ci.org/marcusps/ExpmV.jl.svg?branch=master)](https://travis-ci.org/marcusps/ExpmV.jl)

This is a Julia translation of the MATLAB implementation of Al-Mohy and Higham's
function for computing `expm(t*A)*v` when `A` is sparse, without explicitly computing `expm(A)`.

The original code can be found at [Matlabcentral File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector/all_files), and the theory is explained in the following article:

*Computing the Action of the Matrix Exponential, with an Application to Exponential Integrators*, Awad H. Al-Mohy and Nicholas J. Higham, SIAM Journal on Scientific Computing 2011 33:2, 488-511. ([preprint](http://eprints.ma.man.ac.uk/1426/))

## Installation

Install into Julia using the package manager:

```julia
Pkg.clone("git@github.com:marcusps/ExpmV.jl.git", "ExpmV")
```

## Usage

```julia
expmv(t,A,v)
```

## Benchmarks

Although both the `ExpmV.jl` and `Expokit.jl` implementations are in the early stages of development (`ExpmV.jl` is a direct translation of MATLAB code, and `Expokit.jl` is not fully optimized), here are some crude benchmarks (using `Benchmark.jl`) that indicate large gains over the dense `expm`. Source can be found in `test\benchmark.jl`.

**Benchmark 1**: Complex matrix with density of 0.0191, dimension 100, 100 trials

| Row | Function             | Average     | Relative | Replications |
|-----|----------------------|-------------|----------|--------------|
| 1   | `Expmv.jl`           | 0.0213879  | 5.24362  | 100          |
| 2   | `Expokit.jl`         | 0.00407885 | 1.0      | 100          |
| 3   | Julia's dense `expm` | 0.0686303  | 16.8259  | 100          |

**Benchmark 2**: Complex matrix with density of 0.1913, dimension 100, 100 trials

| Row | Function             | Average    | Relative | Replications |
|-----|----------------------|------------|----------|--------------|
| 1   | `Expmv.jl`           | 0.00602035 | 1.74028  | 100          |
| 2   | `Expokit.jl`         | 0.00345941 | 1.0      | 100          |
| 3   | Julia's dense `expm` | 0.0275519  | 7.96434  | 100          |

**Benchmark 3**: Complex matrix with density of 0.3651, dimension 100, 100 trials

| Row | Function             | Average    | Relative | Replications |
|-----|----------------------|------------|----------|--------------|
| 1   | `Expmv.jl`           | 0.0105414  | 1.97621  | 100          |
| 2   | `Expokit.jl`         | 0.00533413 | 1.0      | 100          |
| 3   | Julia's dense `expm` | 0.0354484  | 6.64558  | 100          |

Clearly the current `ExpmV.jl` implementation needs to be looked at more carefully (See [Issue #1](https://github.com/marcusps/ExpmV.jl/issues/1) in particular, which appears to cause a bottle neck due to my incomplete implementation)

## License

Released under the [BSD 2-clause license](https://tldrlegal.com/license/bsd-2-clause-license-(freebsd)) used in Al-Mohy and  Higham's original code.


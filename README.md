# ExpmV

This is a Julia translation of the MATLAB implementation of Al-Mohy and Higham's
function for computing `expm(A)*v` when `A` is sparse, without explicitly computing `expm(A)`.

The original code can be found at [Matlabcentral File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector/all_files), and the theory is explained in the following article:

*Computing the Action of the Matrix Exponential, with an Application to Exponential Integrators*, Awad H. Al-Mohy and Nicholas J. Higham, SIAM Journal on Scientific Computing 2011 33:2, 488-511. ([preprint](http://eprints.ma.man.ac.uk/1426/))

## Usage

```julia
expmv(t,A,v)
```

## Benchmarks

Although both the `ExpmV.jl` and `Expokit.jl` implementations are in the early stages of development (`ExpmV.jl` is a direct translation of MATLAB code, and `Expokit.jl` is not fully optimized), here are some crude benchmarks (using `Benchmark.jl`) that indicate large gains over the dense `expm`. Source can be found in `test\benchmark.jl`.

**Benchmark 1**: Complex matrix with density of 0.0218, dimension 100, 100 trials

| Row | Function          | Average     | Relative | Replications |
|-----|-------------------|-------------|----------|--------------|
| 1   | `Expmv.jl`   | 0.000409226 | 1.0      | 100          |
| 2   | `Expokit.jl` | 0.00303413  | 7.4143   | 100          |
| 3   | Julia's dense `expm`    | 0.0237006   | 57.9157  | 100          |

**Benchmark 2**: Complex matrix with density of 0.1953, dimension 100, 100 trials

| Row | Function          | Average    | Relative | Replications |
|-----|-------------------|------------|----------|--------------|
| 1   | `Expmv.jl`   | 0.00347449 | 1.0      | 100          |
| 2   | `Expokit.jl` | 0.00605936 | 1.74395  | 100          |
| 3   | Julia's dense `expm`    | 0.0275992  | 7.94337  | 100          |

**Benchmark 3**: Complex matrix with density of 0.3591, dimension 100, 100 trials

| Row | Function          | Average    | Relative | Replications |
|-----|-------------------|------------|----------|--------------|
| 1   | `Expmv.jl`   | 0.0969487  | 9.92966  | 100          |
| 2   | `Expokit.jl` | 0.00976355 | 1.0      | 100          |
| 3   | Julia's dense `expm`    | 0.026019   | 2.66491  | 100          |

Clearly the current `ExpmV.jl` implementation needs to be looked at more carefully in high density cases, but at low densities it performs well.

## License

Released under the [BSD 2-clause license](https://tldrlegal.com/license/bsd-2-clause-license-(freebsd)) used in Al-Mohy and  Higham's original code.

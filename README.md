# ExpmV

[![Build Status](https://travis-ci.org/marcusps/ExpmV.jl.svg?branch=master)](https://travis-ci.org/marcusps/ExpmV.jl)

This is a Julia translation of the MATLAB implementation of Al-Mohy and Higham's
function for computing `expm(A)*v` when `A` is sparse, without explicitly computing `expm(A)`.

The original code can be found at [Matlabcentral File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector/all_files), and the theory is explained in the following article:

*Computing the Action of the Matrix Exponential, with an Application to Exponential Integrators*, Awad H. Al-Mohy and Nicholas J. Higham, SIAM Journal on Scientific Computing 2011 33:2, 488-511. ([preprint](http://eprints.ma.man.ac.uk/1426/))

The fast 1-norm estimation that is crucial for the speed of this algorithm is adapted from code available in Julia's `Base` module. This function (`norm1est`) is licensed under the MIT license.

## Installation

Install into Julia using the package manager:

```julia
Pkg.clone("git@github.com:marcusps/ExpmV.jl.git", "ExpmV")
```

## Usage

```julia
julia> A = sprandn(100,100,.1)+im*sprandn(100,100,.1);

julia> v = randn(100)+im*randn(100); v = v/norm(v);

julia> t = rand();

julia> @time u0 = expm(t*full(A))*v;
  0.015917 seconds (123 allocations: 6.264 MB)

julia> @time u1 = expmv(t,A,v);
  0.037030 seconds (31.68 k allocations: 1.205 MB)

julia> norm(u0 - u1)
8.714711608539494
```

## Benchmarks

Crude benchmarks indicate that `ExpmV.jl` is faster than `Expokit.jl`, and significantly faster than the dense matrix multiplication and multiplication available in Julia. Source can be found in `test\benchmark.jl`.

**Benchmark 1**: Complex matrix with density of ~2%, dimension 100, 100 trials

Implementation | abs lower quartile | abs median | abs upper quartile | rel lower quartile | rel median | rel upper quartile |
---------------|--------------------|------------|--------------------|--------------------|------------|--------------------|
`ExpmV.jl`   | 0.00230458962 | 0.00249493934 | 0.00263433469 | 1 | 1 | 1 |
`Expokit.jl` | 0.009129502759999998 | 0.009715173340000001 | .009802634939999999 | 3.641626053305747 | 3.910858019215716 | 3.9713409243429734 | 
`Base.expm`  | 0.026683810000000002 | 0.027274634689999994 | 0.028157486649999997 | 10.151219683606591 | 10.83995947464359 | 11.97934793700928 |

**Benchmark 2**: Complex matrix with density of ~20%, dimension 100, 100 trials

Implementation | abs lower quartile | abs median | abs upper quartile | rel lower quartile | rel median | rel upper quartile |
---------------|--------------------|------------|--------------------|--------------------|------------|--------------------|
`ExpmV.jl`   | .00230458962|.00249493934|.00263433469|1|1|1|
`Expokit.jl` | .009129502759999998|.009715173340000001|.009802634939999999|3.641626053305747|3.910858019215716|3.9713409243429734|
`Base.expm`  | .026683810000000002|.027274634689999994|.028157486649999997|10.151219683606591|10.83995947464359|11.97934793700928|

**Benchmark 3**: Complex matrix with density of ~40%, dimension 100, 100 trials

Implementation | abs lower quartile | abs median | abs upper quartile | rel lower quartile | rel median | rel upper quartile |
---------------|--------------------|------------|--------------------|--------------------|------------|--------------------|
`ExpmV.jl`   | .0026873555699999997|.0027783622599999998|.0028635493899999997|1|1|1|
`Expokit.jl` | .00976932818|.010370690230000002|.01076657426|3.6008344786481534|3.7028479574276485|3.7624224879752255|
`Base.expm`  | .04804835879999999|.04856447126999999|.049459385209999984|16.963073178252024|17.27205592567062|18.123495132802596|


## License

Released under the [BSD 2-clause license](https://tldrlegal.com/license/bsd-2-clause-license-(freebsd)) used in Al-Mohy and  Higham's original code, although one key component (`norm1est`) is licensed under MIT license, since it is adapted from the implementation in the `Base` module of julia.


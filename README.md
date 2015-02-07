# ExpmV

**Please use the [Github repo](https://github.com/marcusps/ExpmV.jl)** 

This is a Julia translation of the MATLAB implementation of Al-Mohy and Higham's
function for computing `expm(A)*v` when `A` is sparse, without explicitly computing `expm(A)`.

The original code can be found at [Matlabcentral File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector/all_files), and the theory is explained in the following article:

*Computing the Action of the Matrix Exponential, with an Application to Exponential Integrators*, Awad H. Al-Mohy and Nicholas J. Higham, SIAM Journal on Scientific Computing 2011 33:2, 488-511. ([preprint](http://eprints.ma.man.ac.uk/1426/))

## Usage

```julia
expmv(t,A,v)
```

## License

Released under the [BSD 2-clause license](https://tldrlegal.com/license/bsd-2-clause-license-(freebsd)) used in Al-Mohy and  Higham's original code.
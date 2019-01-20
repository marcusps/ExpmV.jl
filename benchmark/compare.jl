using BenchmarkTools, SparseArrays, LinearAlgebra
import ExpmV, Expokit

#BenchmarkTools.DEFAULT_PARAMETERS.seconds = 120

const p = 1e-3
SUITE = BenchmarkGroup()

SUITE["Expm"] = BenchmarkGroup()
SUITE["ExpmV"] = BenchmarkGroup()
SUITE["Expokit"] = BenchmarkGroup()

dimensions = [d for d in 2 .^ (5:10)]

for d in dimensions
    #SUITE[d] = BenchmarkGroup()
    SUITE["Expm"][d] = @benchmarkable exp(full_r)*rv setup=((full_r,rv) = (Matrix(randn()*sprandn($d,$d,p/2)+1im*sprandn($d,$d,p/2)), normalize(randn($d)+1im*randn($d))))
    SUITE["ExpmV"][d] = @benchmarkable ExpmV.expmv(rt,r,rv) setup=((rt,r,rv)=(rand(), sprandn($d,$d,p/2)+1im*sprandn($d,$d,p/2), normalize(randn($d)+1im*randn($d))))
    SUITE["Expokit"][d] = @benchmarkable Expokit.expmv(rt,r,rv) setup=((rt,r,rv)=(rand(), sprandn($d,$d,p/2)+1im*sprandn($d,$d,p/2), normalize(randn($d)+1im*randn($d))))
end

results = run(SUITE, verbose = true)

using BenchmarkTools, ExpmV

const SUITE = BenchmarkGroup()

dimensions = 2 .^ (5:10)
density = 10. .^ (-4 : -1)

for p in density
    SUITE[p] = BenchmarkGroup(["density"])
    for d in dimensions
        SUITE[p][d] = @benchmarkable ExpmV.expmv(rt,r,rv) setup=((rt,r,rv)=(rand(), sprandn($d,$d,$p/2)+1im*sprandn($d,$d,$p/2), normalize(randn($d)+1im*randn($d))))
    end
end

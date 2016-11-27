using BenchmarkTools, Expokit, ExpmV

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 120

const d = 1_000
const p = 1e-4

b1 = @benchmarkable expm(full_r)*rv setup=((full_r,rv) = (full(randn()*sprandn(d,d,p/2)+1im*sprandn(d,d,p/2)), rv=normalize(randn(d)+1im*randn(d))))
b2 = @benchmarkable Expokit.expmv(rt,r,rv) setup=((rt,r,rv)=(rand(), sprandn(d,d,p/2)+1im*sprandn(d,d,p/2), normalize(randn(d)+1im*randn(d))))
b3 = @benchmarkable ExpmV.expmv(rt,r,rv) setup=((rt,r,rv)=(rand(), sprandn(d,d,p/2)+1im*sprandn(d,d,p/2), normalize(randn(d)+1im*randn(d))))

t1 = run(b1)
t2 = run(b2)
t3 = run(b3)

println("Full expm")
println(t1)

println("Expokit")
println(t2)

println("ExpmV")
println(t3)

println("Full expm vs. Expokit")
println(judge(median(t2),median(t1)))

println("Full expm vs. ExpmV")
println(judge(median(t3),median(t1)))

println("Expokit vs. ExpmV")
println(judge(median(t3),median(t2)))

JLD.save("optimize-bench.jld","t1",t1,"t2",t2,"t3",t3)
using BenchmarkTools, Expokit, ExpmV

const d = 100
const p = 0.01

t1 = @benchmark expm(full(randn()*(sprandn(d,d,p/2)+1im*sprandn(d,d,p/2))))*normalize(randn(d)+1im*randn(d))
t2 = @benchmark Expokit.expmv(randn(),sprandn(d,d,p/2)+1im*sprandn(d,d,p/2),normalize(randn(d)+1im*randn(d)))
t3 = @benchmark ExpmV.expmv(randn(),sprandn(d,d,p/2)+1im*sprandn(d,d,p/2),normalize(randn(d)+1im*randn(d)))

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

using BenchmarkTools, Expokit, ExpmV

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 120

const p = 1e-3

bench = Dict{String,Any}[]
for d in 2.^(5:10)
    println("d = $d ...")
    b1 = @benchmarkable expm(full_r)*rv setup=((full_r,rv) = (full(randn()*sprandn($d,$d,p/2)+1im*sprandn($d,$d,p/2)), rv=normalize(randn($d)+1im*randn($d))))
    b2 = @benchmarkable Expokit.expmv(rt,r,rv) setup=((rt,r,rv)=(rand(), sprandn($d,$d,p/2)+1im*sprandn($d,$d,p/2), normalize(randn($d)+1im*randn($d))))
    b3 = @benchmarkable ExpmV.expmv(rt,r,rv) setup=((rt,r,rv)=(rand(), sprandn($d,$d,p/2)+1im*sprandn($d,$d,p/2), normalize(randn($d)+1im*randn($d))))

    t1 = run(b1)
    t2 = run(b2)
    t3 = run(b3)

    push!(bench,Dict("d"=>d, "p"=>p,"expm"=>t1,"expokit"=>t2,"expmv"=>t3))
end

JLD.save("master-bench.jld","bench",bench)

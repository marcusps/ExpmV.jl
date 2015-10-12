using Benchmark, Expokit, ExpmV, DataFrames, DataFramesMeta

const d = 100
const N = 100

if length(ARGS) < 1
    error("Please provide the target matrix density from the command line.")
end

p=0.01

function setup(d,p)
  r = sprandn(d,d,p)+1im*sprandn(d,d,p);

  rv = randn(d)+1im*randn(d); 
  rv = rv/norm(rv,2);

  rt = randn();

  full_r = full(rt*r);  
  return rt,r,rv,full_r
end

function expokitf(rt,r,rv,full_r)
  Expokit.expmv(rt,r,rv)
end

function expmvf(rt,t,rv,full_r)
  ExpmV.expmv(rt,r,rv)
end

function expmf(rt,t,rv,full_r)
  expm(full_r)*rv
end


begin 
    println("Setup ...")
    p = parse(Float64,ARGS[1])
    rt,r,rv,full_r = setup(d,p)
    
    expmv_example() = expokitf(rt,r,rv,full_r)
    expokit_example() = expmvf(rt,r,rv,full_r)
    expm_example() = expmf(rt,r,rv,full_r)

    println("Warming up ...")
    expmv_example()
    expokit_example()
    expm_example()
    
    println("Benchmarking...")

    println("density of $(nnz(r)/prod(size(r))), dimension $d, $N trials")
    t = compare([expmv_example, expokit_example, expm_example],N)

    for i in 1:20
        rt,r,rv,full_r = setup(d,p)
        println("density of $(nnz(r)/prod(size(r))), dimension $d, $N trials")
        append!(t,compare([expmv_example, expokit_example, expm_example],N))
    end

    expmv_q   = quantile(@where(t, :Function .== "expmv_example")[:Average],[.25,.5,.75])
    expokit_q = quantile(@where(t, :Function .== "expokit_example")[:Average],[.25,.5,.75])
    expm_q    = quantile(@where(t, :Function .== "expm_example")[:Average],[.25,.5,.75])

    expmv_rq   = quantile(@where(t, :Function .== "expmv_example")[:Relative],[.25,.5,.75])
    expokit_rq = quantile(@where(t, :Function .== "expokit_example")[:Relative],[.25,.5,.75])
    expm_rq    = quantile(@where(t, :Function .== "expm_example")[:Relative],[.25,.5,.75])

    res = vcat( hcat( expmv_q',   expmv_rq' ),
                hcat( expm_q',    expm_rq' ),
                hcat( expokit_q', expokit_rq' ) ) 

    writecsv("benchmark-p-$(round(Int,p*100)).csv", res)
end


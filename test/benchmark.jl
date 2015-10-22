using Benchmark, Expokit, ExpmV

const d = 100
const N = 100

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

println("Setup ...")
const rt,r,rv,full_r = setup(d,parsefloat(ARGS[1]))

expmv_example() = expmvf(rt,r,rv,full_r)
expokit_example() = expokitf(rt,r,rv,full_r)
expm_example() = expmf(rt,r,rv,full_r)

println("Warming up ...")
expmv_example()
expokit_example()
expm_example()

println("Benchmarking...")
println("density of $(nnz(r)/prod(size(r))), dimension $d, $N trials")
c=compare([expmv_example, expokit_example, expm_example],N)
println(c)



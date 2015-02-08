using Benchmark, Expokit, ExpmV

d = 100
N = 100

function setup(d,p)
  r = sprandn(d,d,p)+1im*sprandn(d,d,p);

  rv = randn(d)+1im*randn(d); 
  rv = rv/norm(rv,2);

  rt = randn();

  full_r = full(rt*r);  
  return rt,r,rv,full_r
end

function expokit_example()
  Expokit.expmv(rt,r,rv)
end

function expmv_example()
  ExpmV.expmv(rt,r,rv)
end

function expm_example()
  expm(full_r)*rv
end

rt,r,rv,full_r = setup(d,.01)
expmv_example()
expokit_example()
expm_example()

rt,r,rv,full_r = setup(d,.01)
println("density of $(nnz(r)/prod(size(r))), dimension $d, $N trials")
c10=compare([expmv_example, expokit_example, expm_example],N)
println(c10)

rt,r,rv,full_r = setup(d,.1)
println("density of $(nnz(r)/prod(size(r))), dimension $d, $N trials")
c20=println(compare([expmv_example, expokit_example, expm_example],N))
println(c20)

rt,r,rv,full_r = setup(d,.2)
println("density of $(nnz(r)/prod(size(r))), dimension $d, $N trials")
c30=println(compare([expmv_example, expokit_example, expm_example],N))
println(c30)


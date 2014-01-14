using Benchmark, ExpmV

d = 100
N = 1000

r = sprandn(d,d,.1)+1im*sprandn(d,d,.1);

rv = randn(d)+1im*randn(d); 
rv = rv/norm(rv,2);

rt = randn();

full_r = full(rt*r);  

function expmv_example()
  (x,_,_,_,_,_) = expmv(rt,r,rv)
  x
end

function expm_example()
  expm(full_r)*rv
end

println(compare([expmv_example, expm_example],N))


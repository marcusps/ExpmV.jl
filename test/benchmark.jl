using Benchmark, ExpmV

d = 100
N = 1000

r = sprandn(d,d,.1)+1im*sprandn(d,d,.1);

rv = randn(d)+1im*randn(d); 
rv = rv/norm(rv,2);

rt = randn();

full_r = full(rt*r);  

function f1()
  (x,_,_,_,_,_) = expmv(rt,r,rv)
  x
end

function f2()
  expm(rt*full_r)*rv
end

compare([f1,f2],N)


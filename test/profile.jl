using ExpmV

const d = 100
const N = 100

p=.02

function setup(d,p)
  r = sprandn(d,d,p)+1im*sprandn(d,d,p);

  rv = randn(d)+1im*randn(d); 
  rv = rv/norm(rv,2);

  rt = randn();

  full_r = full(rt*r);  
  return rt,r,rv,full_r
end

function expmvf(rt,t,rv,full_r)
    for i = 1:1000
        expmv(rt,r,rv)
    end
end

println("Setup ...")
const rt,r,rv,full_r = setup(d,p)

expmv_example() = expmvf(rt,r,rv,full_r)

println("Warming up ...")
expmv_example()

@profile expmv_example()



using ExpmV
using Base.Test

for i = 1:100, d in 10:10:60, herm in [true, false]

    d = 60

    r = sprandn(d,d,.1)+1im*sprandn(d,d,.1)
    if herm 
        r = (r-r')/2
    end

    rv = randn(d)+1im*randn(d)
    rv = rv / norm(rv,2)
    
    rt = randn()

    x = expmv(rt,r,rv)

    @test_approx_eq_eps norm(x-expm(full(rt*r))*rv,2) 0.0 1e-9
end

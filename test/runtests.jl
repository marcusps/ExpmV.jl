using ExpmV
using Base.Test

for i = 1:100, d in 10:10:60, herm in [true, false]
    r = sprandn(d,d,.1)+1im*sprandn(d,d,.1)
    if herm
        r = (r-r')/2
    end

    rv = randn(d)+1im*randn(d)
    rv = rv / norm(rv,2)

    rt = randn()

    x = expmv(rt,r,rv)
    @test norm(x-expm(full(rt*r))*rv,2) ≈ 0.0 atol=1.0e-9

    # Test the StepRangeLen version against the normal version

    t = linspace(0, rt, 51)
    x = expmv(t,r,rv)
    y = hcat([expmv(ti,r,rv) for ti in t]...)
    @test x ≈ y atol=1.0e-10
end

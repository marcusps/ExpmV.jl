using ExpmV
using Test
using LinearAlgebra
using SparseArrays

@testset "Hermitian: $herm"  for herm in [true, false]
    @testset "Size: $d" for d in 10:10:60
        for i = 1:20
            r = sprandn(d,d,.1)+1im*sprandn(d,d,.1)
            if herm
                r = (r-r')/2
            end

            rv = randn(d)+1im*randn(d)
            rv = rv / norm(rv,2)

            rt = randn()

            x = expmv(rt,r,rv)
            @testset "Against expm" begin
                @test norm(x-exp(Matrix(rt*r))*rv,2) ≈ 0.0 atol=1.0e-9
            end

            # Test the StepRangeLen version against the normal version
            @testset "Timespan $nt timesteps" for nt in [5 11 51]
                t = range(0, stop=rt, length=nt)
                x = expmv(t,r,rv)
                y = hcat([expmv(ti,r,rv) for ti in t]...)
                @test x ≈ y atol=1.0e-9
            end

            @testset "Second dim: $d2" for d2 in 2:4
                rv = randn(d,d2)+1im*randn(d,d2)
                rv = rv*diagm(0 => [1.0/norm(rv[:,j],2) for j in 1:d2])

                x = expmv(rt,r,rv)
                @test norm(x-exp(Matrix(rt*r))*rv,2) ≈ 0.0 atol=1.0e-9
            end
        end
    end
end

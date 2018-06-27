using ExpmV

if VERSION < v"0.7-"
    using Base.Test
else
    using Test
    using LinearAlgebra
    using SparseArrays
end

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
                if VERSION < v"0.7-"
                    t = linspace(0, rt, nt)
                else
                    t = range(0, stop=rt, length=nt)
                end
                x = expmv(t,r,rv)
                y = hcat([expmv(ti,r,rv) for ti in t]...)
                @test x ≈ y atol=1.0e-10
            end
        end
    end
end

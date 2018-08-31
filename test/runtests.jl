using ExpmV
using Test
using LinearAlgebra
using SparseArrays

@testset "Real matrices" begin
    @testset "Positive: $positive"  for positive in [true, false]
        @testset "Size: $d" for d in 10:20:70
            for i = 1:10
                r = sprandn(d, d, .1)

                if positive
                    r += abs.(r)
                end

                rv = randn(d)
                rv = rv / norm(rv,2)

                rt = randn()

                x = expmv(rt,r,rv)

                @testset "Against expm" begin
                    @test x ≈ exp(Matrix(rt*r))*rv
                end

                # Test the StepRangeLen version against the normal version
                @testset "Timespan $nt timesteps" for nt in [5 11 51]
                    t = range(0, stop=rt, length=nt)
                    x = expmv(t,r,rv)
                    y = hcat([expmv(ti,r,rv) for ti in t]...)
                    @test x ≈ y
                end

                @testset "Second dim: $d2" for d2 in 2:4
                    rv = randn(d,d2)+1im*randn(d,d2)
                    rv = rv*diagm(0 => [1.0/norm(rv[:,j],2) for j in 1:d2])

                    x = expmv(rt,r,rv)
                    @test x ≈ exp(Matrix(rt*r))*rv
                end
            end
        end
    end
end

@testset "Complex matrices" begin
    @testset "Hermitian: $herm"  for herm in [true, false]
        @testset "Size: $d" for d in 10:20:70
            for i = 1:10
                r = sprandn(d,d,.1)+1im*sprandn(d,d,.1)
                if herm
                    r = (r-r')/2
                end

                rv = randn(d)+1im*randn(d)
                rv = rv / norm(rv,2)

                rt = randn()

                x = expmv(rt,r,rv)
                @testset "Against expm" begin
                    @test x ≈ exp(Matrix(rt*r))*rv
                end

                # Test the StepRangeLen version against the normal version
                @testset "Timespan $nt timesteps" for nt in [5 11 51]
                    t = range(0, stop=rt, length=nt)
                    x = expmv(t,r,rv)
                    y = hcat([expmv(ti,r,rv) for ti in t]...)
                    @test x ≈ y
                end

                @testset "Second dim: $d2" for d2 in 2:4
                    rv = randn(d,d2)+1im*randn(d,d2)
                    rv = rv*diagm(0 => [1.0/norm(rv[:,j],2) for j in 1:d2])

                    x = expmv(rt,r,rv)
                    @test x ≈ exp(Matrix(rt*r))*rv
                end
            end
        end
    end
end

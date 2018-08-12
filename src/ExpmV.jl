module ExpmV

  using LinearAlgebra
  using SparseArrays
  # package code goes here
  include("expmv_fun.jl")
  include("degree_selector.jl")
  include("normAm.jl")
  include("expmv_tspan.jl")
  include("select_taylor_degree.jl")

end # module

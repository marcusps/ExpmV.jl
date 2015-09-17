VERSION >= v"0.4.0-dev+6521" && __precompile__()
module ExpmV

  # package code goes here
  include("normAm.jl")

  include("select_taylor_degree.jl")

  include("expmv_fun.jl")

end # module

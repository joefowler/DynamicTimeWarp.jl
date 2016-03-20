
module DynamicTimeWarp

using Distances

export dtw,
       dba

include("dtw.jl")
include("dba.jl")
#include("windowed_matrix.jl")
#include("fastdtw.jl")

end # module

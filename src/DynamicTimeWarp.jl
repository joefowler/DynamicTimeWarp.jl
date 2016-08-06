
module DynamicTimeWarp

using Distances

export dtw,
       dba,
       fastdtw

include("dtw.jl")
include("dba.jl")
include("windowed_matrix.jl")
include("fastdtw.jl")

end # module

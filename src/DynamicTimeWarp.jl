module DynamicTimeWarp

# package code goes here

# Dynamic Time Warping with a user-specified distance function

function dtw(seq1::Vector, seq2::Vector, distance::Function=Distance.square)
    const m=length(seq1)
    const n=length(seq2)

    cost11 = distance(seq1[1], seq2[1])
    cost = zeros(typeof(cost11), m, n)

    # Initialize first column and first row
    cost[1,1] = cost11
    for r=2:m
        cost[r,1] = cost[r-1,1] + distance(seq1[r], seq2[1])
    end
    for c=2:n
        cost[1,c] = cost[1,c-1] + distance(seq1[1], seq2[c])
    end

    # Build cost matrix
    for c=2:n
        for r=2:m
            best_neighbor_cost = min(cost[r-1,c], cost[r-1,c-1], cost[r,c-1])
            cost[r,c] = best_neighbor_cost + distance(seq1[r], seq2[c])
        end
    end

    match1, match2 = trackback(cost)
    
    cost[m,n], match1, match2
end

function trackback(D)
    i,j = size(D)
    p,q = [i],[j]
    while i > 1 && j > 1
        tb = indmin([D[i-1,j-1], D[i-1,j], D[i,j-1]])
        tb in [1,2] && (i-=1)
        tb in [1,3] && (j-=1)
        unshift!(p,i)
        unshift!(q,j)
    end
    unshift!(p,1)
    unshift!(q,1)
    p,q
end

using PyPlot

function plotdtw{T<:Real}(seq1::Vector{T}, seq2::Vector{T}, offset::Float64=0.0)
    cost, match1, match2 = dtw(seq1, seq2)
    if offset == 0.0
        offset = 2*(std(seq1) + std(seq2)) + mean(seq1) - mean(seq2)
    end
    clf()
    plot(seq1, "-r", seq2+offset, "-b")
    for i=1:length(match1)
        plot([match1[i],match2[i]]-1, [seq1[match1[i]], seq2[match2[i]]+offset],
             color="gray")
    end
end



# This sub-module contains various functions returning pointwise distances
# between elements from two seequences.

module Distance
square(x,y) = (x-y)^2
absval(x,y) = abs(x-y)

# And now a Poisson-sensitive distance.
# Call this with the 2 sequences, and it returns a closure, a function that
# measures Poisson-distance for elements of those two particular sequences
function poissonclosure(seq1, seq2)
    n1=sum(seq1)
    n2=sum(seq2)
    function pdist(x,y)
        ### WARNING! This does not yet weight properly for Poissonianness
        ### It only accounts for differences in normalization.
        (x/n1 - y/n2)^2
    end
end

end # module Distance


export
    dtw

end # module

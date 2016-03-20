function dtw(seq1::Vector, seq2::Vector, distance::Function=Distance.square)

    # Build the cost matrix
    const m=length(seq2)
    const n=length(seq1)
    cost11 = distance(seq1[1], seq2[1])
    cost = zeros(typeof(cost11), m, n)

    # Initialize first column and first row
    cost[1,1] = cost11
    for r=2:m
        cost[r,1] = cost[r-1,1] + distance(seq1[1], seq2[r])
    end
    for c=2:n
        cost[1,c] = cost[1,c-1] + distance(seq1[c], seq2[1])
    end

    # Complete the cost matrix
    for c=2:n
        for r=2:m
            best_neighbor_cost = min(cost[r-1,c], cost[r-1,c-1], cost[r,c-1])
            cost[r,c] = best_neighbor_cost + distance(seq1[c], seq2[r])
        end
    end

    trackcols, trackrows = trackback(cost)
    cost[end,end], trackcols, trackrows
end


"""
cols,rows = trackback(D::Matrix{AbstractFloat})

Given the cost matrix `D`, computes the optimal track from end to beginning.
Returns `cols` and `rows` which are vectors respectively holding the track.
"""
function trackback(D::Matrix{AbstractFloat})
    r,c = size(D)
    rows,cols = Int[r],Int[c]
    while r > 1 && c > 1
        tb = indmin([D[r-1,c-1], D[r-1,c], D[r,c-1]])
        tb in [1,2] && (r-=1)
        tb in [1,3] && (c-=1)
        push!(rows,r)
        push!(cols,c)
    end
    # Possibly either r>1 or c>1 at this point (but not both). 
    # Add the unfinished part of the track to reach [1,1]
    for r=r-1:-1:1
        push!(rows,r)
        push!(cols,1)
    end
    for c=c-1:-1:1
        push!(rows,1)
        push!(cols,c)
    end
    reverse(cols), reverse(rows)
end

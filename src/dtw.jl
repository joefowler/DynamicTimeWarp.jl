"""
    cost,i1,i2 = dtw(seq1, seq2, [dist=SqEuclidean])

Find a set of indices (`i1`,`i2`) that align two time series (`seq1`,`seq2`) by
dynamic time warping. Also returns the distance (after warping) according to
the SemiMetric `dist`, which defaults to squared Euclidean distance (see
Distances.jl). If `seq1` and `seq2` are matrices, each column is considered
an observation.
"""
function dtw(
        seq1::Vector,
        seq2::Vector,
        dist::SemiMetric=SqEuclidean()
    )

    # Build the cost matrix
    m = length(seq2)
    n = length(seq1)
    cost11 = evaluate(dist, seq1[1], seq2[1])
    cost = zeros(typeof(cost11), m, n)

    # Initialize first column and first row
    cost[1,1] = cost11
    for r=2:m
        cost[r,1] = cost[r-1,1] + evaluate(dist, seq1[1], seq2[r])
    end
    for c=2:n
        cost[1,c] = cost[1,c-1] + evaluate(dist, seq1[c], seq2[1])
    end

    # Complete the cost matrix
    for c=2:n
        for r=2:m
            best_neighbor_cost = min(cost[r-1,c], cost[r-1,c-1], cost[r,c-1])
            cost[r,c] = best_neighbor_cost + evaluate(dist, seq1[c], seq2[r])
        end
    end

    trackcols, trackrows = trackback(cost)
    cost[end,end], trackcols, trackrows
end

# Wrapper for multi-dimensional time series.
function dtw(
        seq1::AbstractMatrix,
        seq2::AbstractMatrix,
        dist::SemiMetric = SqEuclidean()
    )
    if size(seq1,1) != size(seq2,1)
        throw(ArgumentError("Provided time series don't have same number of variables, size(seq1,1) must equal size(seq2,1)"))
    end
    s1 = [view(seq1,:,i) for i = 1:size(seq1,2)]
    s2 = [view(seq2,:,i) for i = 1:size(seq2,2)]
    dtw(s1,s2,dist)
end

"""
    cost,i1,i2 = dtw(seq1,seq2,i2min,i2max,[dist=SqEuclidean])

Do DTW to align `seq1` and `seq2` confined to a window. Vectors `i2min` and
`i2max` specify (inclusive) lower and upper bounds for `seq2` for each index in
`seq1`. Thus, `i2min` and `i2max` are required to be the same length as `seq1`.
"""
function dtw(seq1::Vector, seq2::Vector,
             i2min::Vector, i2max::Vector,
             dist::SemiMetric = SqEuclidean())

    m = length(seq2) # of rows  in cost matrix
    n = length(seq1) # of columns in cost matrix
    n == length(i2min) || throw(ArgumentError("i2min does not match length of seq1."))
    n == length(i2max) || throw(ArgumentError("i2max does not match length of seq1."))
    1 == i2min[1] || throw(ArgumentError("i2min must start at 1."))
    m == i2max[end] || throw(ArgumentError("i2max must end at length(seq2)."))

    # Build the (n x m) cost matrix into a WindowedMatrix, because it's ragged.
    # That type gives efficient storage with convenient [r,c] indexing and returns
    # Inf when accessed outside the window.
    cost = WindowedMatrix(i2min, i2max, Inf)

    # First column first
    cost[1,1] = evaluate(dist, seq1[1], seq2[1])
    for r=2:i2max[1]
        cost[r,1] = cost[r-1,1]  + evaluate(dist, seq1[1], seq2[r])
    end

    # Complete the cost matrix from columns 2 to m.
    for c=2:n
        for r=i2min[c]:i2max[c]
            best_neighbor_cost = min(cost[r-1,c], cost[r-1,c-1], cost[r,c-1])
            cost[r,c] = best_neighbor_cost + evaluate(dist, seq1[c], seq2[r])
        end
    end
    trackcols, trackrows = trackback(cost)
    cost[end,end], trackcols, trackrows
end


"""
    cols,rows = trackback(D::Matrix)

Given the cost matrix `D`, computes the optimal track from end to beginning.
Returns `cols` and `rows` which are vectors respectively holding the track.
"""
function trackback{T<:Number}(D::AbstractMatrix{T})
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

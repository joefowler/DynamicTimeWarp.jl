# The FastDTW approximation to the DTW, described in "FastDTW: Toward Accurate
# Dynamic Time Warping in Linear Time and Space", S Salvador & P Chan, __Intelligent
# Data Analysis__ (2007).

function fastdtw(seq1::Vector, seq2::Vector, radius::Integer, 
                 distance::Function=Distance.square)
    const MinSize = max(radius + 2, 10)
    const N1 = length(seq1)
    const N2 = length(seq2)
    if N1 <= MinSize || N2 <= MinSize
        return (dtw(seq1, seq2, distance))
    end

    # Call recursively on a pair of sequences half this length
    compressed1 = compress(seq1)
    compressed2 = compress(seq2)
    _cost, lowrescol, lowresrow = fastdtw(compressed1, compressed2, radius, distance)

    # Now resample that path to the finer resolution, find the correct
    # window around it, and get the DTW given that window.
    hirescol, hiresrow = expandpath(lowrescol, lowresrow, N1, N2)
    idx2min, idx2max = computewindow(hirescol, hiresrow, radius)
    cost1, newcol, newrow = dtwwindowed(seq1, seq2, idx2min, idx2max, distance)
end



# Given a path through low-res space, generate an approximate path
# through high-res space. It should have dimension Ncol x Nrow

function expandpath(lowrescol, lowresrow, Ncol, Nrow)
    @assert div(Ncol+1,2) == lowrescol[end]
    @assert div(Nrow+1,2) == lowresrow[end]
    const Np = length(lowrescol)
    @assert Np == length(lowresrow)

    hirescol = zeros(eltype(lowrescol), 2*Np)
    hiresrow = zeros(eltype(lowresrow), 2*Np)
    hirescol[1] = hiresrow[1] = c = r = 1
    for i=1:Np-1
        # Select plan according to the next move in lowres path.
        if lowrescol[i+1] == lowrescol[i]  # Next move is up
            r += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            r += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r
            
        elseif lowresrow[i+1] == lowresrow[i] # Next move is sideways
            c += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            c += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r
            
        else  # Next move is diagonal.
            c += 1; r += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            c += 1; r += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r
        end
    end
    hirescol[end] = Ncol
    hiresrow[end] = Nrow
    # When expanding to an odd numbered size, it's possible to repeat
    # the last step.  Fix that:
    if hirescol[end]==hirescol[end-1] && hiresrow[end]==hiresrow[end-1]
        hirescol = hirescol[1:end-1]
        hiresrow = hiresrow[1:end-1]
    end
    hirescol, hiresrow
end


# Given the lists of (col,row) indices for the optimal path, compute a "window"
# around that path of the given radius.
# Returns (rowmin, rowmax), each a vector of length pathcols[end], representing
# for each column, the minimum and maximum row numbers used in that column.

function computewindow(pathcols, pathrows, radius)
    const Np = length(pathcols)
    @assert Np == length(pathrows)
    const Ncol = pathcols[end]
    const Nrow = pathrows[end]

    # Find the min/max row at each column in the path.
    pathmin = zeros(Int, Ncol)
    pathmax = zeros(Int, Ncol)
    for i=1:Np
        c,r = pathcols[i], pathrows[i]
        pathmax[c] = r
        if pathmin[c] == 0
            pathmin[c] = r
        end
    end

    # The window in each column for "radius" r starts at the pathmin
    # of the rth-previous column and ends at the pathmax of the
    # rth-next column, plus (in each case) the radius.
    if radius < Ncol-1 && radius < Nrow-1
        rowmin = vcat(fill(1,radius), pathmin[1:end-radius]-radius)
        rowmax = vcat(pathmax[radius+1:end]+radius, fill(Nrow,radius))

        # Window values must be in the range [1:Nrow].
        for c=1:Ncol
            if rowmin[c]<1; rowmin[c]=1; end
            if rowmax[c]>Nrow; rowmax[c]=Nrow; end
        end
    else
        rowmin = fill(1,Ncol)
        rowmax = fill(Nrow,Ncol)
    end
    rowmin, rowmax
end



# Do DTW in a subset of the full space, the subset specified by
# a "window". Arguments [idx2min,idx2max] give the inclusive lowest
# and highest index in the seq2 direction, one element for each index along 
# the seq1 direction. Thus seq1, idx2min, and idx2max should all be of
# equal length

function dtwwindowed(seq1::Vector, seq2::Vector,
                     idx2min::Vector, idx2max::Vector,
                     distance::Function=Distance.square)

    const m=length(seq2) # of rows  in cost matrix
    const n=length(seq1) # of columns in cost matrix
    @assert n==length(idx2min)
    @assert n==length(idx2max)
    @assert 1==minimum(idx2min)
    @assert m==maximum(idx2max)

    # Build the (n x m) cost matrix into a WindowedMatrix, because it's ragged.
    # That type gives efficient storage with convenient [r,c] indexing and returns
    # Inf when accessed outside the window.
    cost = WindowedMatrix(idx2min, idx2max, Inf)

    # First column first
    cost[1,1] = distance(seq1[1], seq2[1])
    for r=2:idx2max[1]
        cost[r,1] = cost[r-1,1]  + distance(seq1[1], seq2[r])
    end

    # Complete the cost matrix from columns 2 to m.
    for c=2:n
        for r=idx2min[c]:idx2max[c]
            best_neighbor_cost = min(cost[r-1,c], cost[r-1,c-1], cost[r,c-1])
            cost[r,c] = best_neighbor_cost + distance(seq1[c], seq2[r])
        end
    end
    trackcols, trackrows = trackback(cost)
    cost[end,end], trackcols, trackrows
end

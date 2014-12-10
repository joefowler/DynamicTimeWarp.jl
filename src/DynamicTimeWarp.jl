module DynamicTimeWarp

# package code goes here

# Dynamic Time Warping with a user-specified distance function

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


# Compute the optimal track backwards through the cost matrix from end to beginning.
# Return (columns, rows) of the optimal track.

function trackback(D)
    r,c = size(D)
    rows,cols = [r],[c]
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


# Do DTW in a subset of the full space, the subset specified by
# a "window". Arguments [idx2min,idx2max] give the inclusive lowest
# and highest index in the seq2 direction, one element for each index along 
# the seq1 direction. Thus seq1, idx2min, and idx2max should all be of
# equal length

function dtwwindowed(seq1::Vector, seq2::Vector,
                     idx2min::Vector, idx2max::Vector,
                     distance::Function=Distance.square)
    const m=length(seq2)
    const n=length(seq1)
    @assert n==length(idx2min)
    @assert n==length(idx2max)

    # Build the (n x m) cost matrix into a 1D vector, because it's ragged
    # The order of elements will be all elements with i (seq1 index) = 1,
    # from j (seq2 index) = idx2min[1]:idx2max[1]. After that, all of i=2, ...
    
    nseq2 = idx2max + 1 - idx2min               # Window height per column
    const N = sum(nseq2)                        # Total # of path steps to test
    idxcol = 1+vcat([0],cumsum(nseq2[1:end-1])) # index of 1st element in each column
    
    cost11 = distance(seq1[1], seq2[1])
    cost = zeros(typeof(cost11), N)

    # First column first
    cost[1] = cost11
    for r=2:idx2max[1]
        cost[r] = cost[r-1]  + distance(seq1[1], seq2[r])
    end

    # Complete the cost matrix from columns 2 to m.
    for c=2:m

        # Handle the lowest element in the column separately
        r = idx2min[c]
        i = idxcol[c]
        idxleft = idxcol[c-1] + idx2min[c]-idx2min[c-1]
        if r<=idx2min[c-1] # Left but no diagonal neighbor
            cost[i] = cost[idxleft] + distance(seq1[c], seq2[r])
        else  # Left and diagonal neighbors
            cost[i] = min(cost[idxleft],cost[idxleft-1]) + distance(seq1[c], seq2[r])
        end

        # Now all other elements in the column
        for r=idx2min[c]+1:idx2max[c]
            i += 1
            idxleft += 1
            if r <= idx2max[c-1] # Left, diag, and lower neighbors
                cost[i] = min(cost[idxleft],cost[idxleft-1],cost[i-1]
                              ) + distance(seq1[c], seq2[r])

            elseif r == idx2max[c-1]+1  # Diag and lower neighbors
                cost[i] = min(cost[idxleft-1], cost[i-1]) + distance(seq1[c], seq2[r])

            else # Only a lower neighbor
                cost[i] = cost[i-1] + distance(seq1[c], seq2[r])
            end
        end
    end

    trackcols, trackrows = trackbackwindowed(cost, idx2min, idx2max)
    cost[end], trackcols, trackrows
end


# Compute the optimal track backwards through the cost matrix from end to beginning.
# Return (columns, rows) of the optimal track.

function trackbackwindowed(D::Vector, idx2min::Vector, idx2max::Vector)
    r = idx2max[end]
    c = length(idx2max)  # r rows, c columns

    nseq2 = idx2max + 1 - idx2min               # Window height per column
    idxcol = 1+vcat([0],cumsum(nseq2[1:end-1])) # index of 1st element in each column

    rows,cols = [r],[c]
    while r > 1 && c > 1
        idx = idxcol[c] + r - idx2min[c]
        Ddown = Inf
        if r > idx2min[c]
            Ddown = D[idx-1]
        end
        Dleft = Inf
        if r >= idx2min[c-1] && r <= idx2max[c-1]
            Dleft = D[idxcol[c-1]+(r-idx2min[c-1])]
        end
        Ddiag = Inf
        if r > idx2min[c-1] && r-1 <= idx2max[c-1]
            Ddiag = D[idxcol[c-1]+(r-1-idx2min[c-1])]
        end

        tb = indmin([Ddiag, Ddown, Dleft])
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






using PyPlot

function plotdtw(seq1::Vector, seq2::Vector, offset=0.0)
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



#function dtwbaryavg_iteration(dbavg::Vector, sequences::Array{Array{1},1})
function dtwbaryavg_iteration(dbavg::Vector, sequences::Array)
    const Nseq = length(sequences)
    count = zeros(Int, length(dbavg))
    sumcoords = zeros(Float64, length(dbavg))
    
    for i=1:Nseq
        cost, match1, match2 = dtw(dbavg, sequences[i])
        for j=1:length(match2)
            count[match1[j]] += 1
            sumcoords[match1[j]] += sequences[i][match2[j]]
        end
        println("Compared $i to the standard")
    end
    sumcoords = sumcoords ./ count
end


#function dtwbaryavg{T<:Real}(sequences::Array{Array{T,N},1})
function dtwbaryavg(sequences)
    dbavg = copy(sequences[1])
    clf()
    for i = 1:5
        dbavg = dtwbaryavg_iteration(dbavg, sequences)
        plot(dbavg)
    end
    dbavg
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
dtw,
dtwwindowed,
dtwbaryavg

end # module

using PyPlot

function dtwspeedtest(N::Integer, functype=1)
    @assert N>10
    if functype == 1
        stretch = div(N,10)
        stretch = max(stretch, 10)
        t = linspace(0,1,N-stretch)
        x = vcat(fill(0.0, stretch), sin((t.^1.3).*5))
        y = vcat(sin((t.^0.7).*5), fill(sin(5), stretch))
    else
        # Random path
        rows=Array(Int, N)
        cols=Array(Int, N)
        r=c=1
        const PSTRETCH=0.8
        for i=1:N
            rnum = rand()
            if rnum > PSTRETCH # go diag
                r+=1; c+=1
            elseif rnum > 0.5*PSTRETCH
                r+=1
            else
                c+=1
            end
            rows[i] = r
            cols[i] = c
        end
        t = linspace(0, N, max(r[end],c[end]))
        f = sin((N-t).^0.7)+ (t.*(3./N))
        x = f[rows]
        y = f[cols]
    end

    clf()
    plt.subplot(211)
    dstd = std(x)+std(y)
    plot(x-mean(x), "r")
    plot(y-mean(y)+dstd, "b")

    plt.subplot(224)
    radii=[40,20,10,5,2]
    colors=["red","orange","gold","green","cyan"]
    for (rad,color) in zip(radii, colors)
        print("Testing size $N radius $rad: ")
        @time _,r1,c1=fastdtw(x,y,rad)
        plot(r1,c1,color=color)
    end
    
    if N<=2500
        @time _,r,c=dtw(x,y)
        plot(r,c,"k")
    end
    nothing
end


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
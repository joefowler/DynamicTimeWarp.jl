using DynamicTimeWarp
using Plots; gr()

function plotdtw(seq1::Vector, seq2::Vector, offset=0.0)
    cost, i1, i2 = dtw(seq1, seq2)
    if offset == 0.0
        offset = 2*(std(seq1) + std(seq2)) + mean(seq1) - mean(seq2)
    end
    plot()
    for (a,b) in zip(i1,i2)
        plot!([a,b], [seq1[a], seq2[b]+offset],
             line=(:grey), label=())
    end
    plot!(legend=(:none))
    plot!(seq1, line=(:red,2))
    plot!(seq2+offset, line=(:blue,2))
end


seq1 = exp(-(linspace(-3,3,100)).^2)
seq2 = vcat(zeros(30),seq1)
plotdtw(seq1,seq2,1.0)

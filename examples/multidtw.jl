using DynamicTimeWarp
using Plots; gr()

function plotdtw(seq1, seq2, offset=0.0)
    cost, i1, i2 = dtw(seq1, seq2)

    # transpose signals for plotting
    s1 = transpose(seq1)
    s2 = transpose(seq2)
    n = size(s1,2) # number of variables

    plot(layout=n)
    for (a,b) in zip(i1,i2)
        plot!([a,b], vcat(s1[a,:], s2[b,:]+offset),
              line=(:grey), label=())
    end
    plot!(s1, line=(:red,2))
    plot!(s2+offset, line=(:blue,2))
    plot!(legend=(:none))
end

f(x) = exp(-x.^2)
g(x) = x.*exp(-x.^2)

t = linspace(-3,3,100)
x1 = f(t)
x2 = vcat(zeros(40),x1)

y1 = g(t)
y2 = vcat(zeros(20),y1,zeros(20))

seq1,seq2 = transpose(hcat(x1,y1)),transpose(hcat(x2,y2))
plotdtw(seq1,seq2,1.0)

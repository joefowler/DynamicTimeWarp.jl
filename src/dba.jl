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

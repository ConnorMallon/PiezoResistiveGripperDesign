using DrWatson
using DataFrames

cols_to_filter = ["Js"]
results  = DrWatson.collect_results(datadir("results/ip_36"))
results  = DrWatson.collect_results(datadir("results/TOIP_72"))

@show names(results)
J3 = results.J3
J1 = results.J1
J2 = results.J2

α1 = results.α1
α2 = results.α2
α3 = results.α3



α1j1α2j2 = - α1 .* J1 + α2 .* J2


using Plots
scatter(α1j1α2j2, -J3,group=α3, xlabel="α1 * J1 + α2 * J2", ylabel="J3", title="J3 vs α1 * J1 + α2 * J2", legend=true)


# Get the sorted indices based on the last value of Js
sorted_indices = sortperm(1:nrow(results), by = i -> last(results.Js[i]))

# Extract top 10 rows
top10 = results[sorted_indices[1:10], :]

bottom10 = results[sorted_indices[end-9:end], :]

last.(top10[:,"Js"])
last.(bottom10[:,"Js"])

top10[:,"Sϕ"]
bottom10[:,"Sϕ"]

top1 = results[sorted_indices[1], :]
top1["Sϕ"]

bottom1 = results[sorted_indices[end], :]
bottom1[:,"Sϕ"]

desired_Sϕ =  [1, 2, 5, 6, 7, 8]
filtered = filter(row -> row.Sϕ == desired_Sϕ, results)
Js = filtered.Js[1]
using Plots
plot(1:length(Js),Js)





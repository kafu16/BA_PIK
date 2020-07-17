include("code_BA.jl")
using StatsPlots, StatsBase, LaTeXStrings

### analysis data final_runs
## decreasing energy for increasing number of steps
# average
Data_av_e_std = collect_av_e_std(SData)
Plots.plot(Data_av_e_std[1], title = L"\textrm{Average } G_{av} \textrm{ depending on number of steps}",
        label = L"\textrm{Average } G_{av}",
        xaxis = L"\textrm{Number of steps}", yaxis = L"\textrm{Average } G_{av}",
        xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, titlefontsize = 10,
        ribbon=Data_av_e_std[2], fillalpha=.2,
        framestyle = :box, grid = true)
Plots.savefig("decreasing_Gav_av.pdf")

# single run
en = Data[8][8]
Plots.plot(en, title = L"\textrm{Decreasing } G_{av} \textrm{ for one sample grid}",
        label = L"G_{av}",
        xaxis = L"\textrm{Number of steps}", yaxis = L"G_{av}",
        xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, titlefontsize = 10,
        framestyle = :box, grid = true)
Plots.savefig("decreasing_Gav_sample.pdf")
Data[6]

## results SA: calculating Gav_av and STD_Gav (averaged over all runs)
# random grids
N_runs = 10
Gav_av = mean(collect_data2(Data, 2, 2))
STD_Gav = 1 / sqrt(N_runs) * std(collect_data2(Data, 2, 2); corrected=true) # corrected=true is default of std(), so it could be omitted
# minimized G_av
Gav_av = mean(collect_data2(Data, 5, 2))
STD_Gav = 1 / sqrt(N_runs) * std(collect_data2(Data, 5, 2); corrected=true)
# maximized G_av
Gav_av = mean(collect_data2(Data_m, 5, 2))
STD_Gav = 1 / sqrt(N_runs) * std(collect_data2(Data_m, 5, 2); corrected=true)

## example grids
# random grid
g = visualize_data(SData, 1, 18)
draw(PDF("random_grid_sample.pdf", 16cm, 16cm), g)
Data[18] #G_av = 85.700
# grid with minimized G_av
g = visualize_data(SData, 4, 9)
draw(PDF("minimized_Gav_grid_sample.pdf", 16cm, 16cm), g)
Data[18] #G_av = 23.006

## flow distributions
# random vs minimized G_av
nbins = 11
random = weights_mean_err(SData, 1, nbins)
minimized = weights_mean_err(SData, 4, nbins)

weights1 = convert(Vector{Float64}, random[1])
weights2 = convert(Vector{Float64}, minimized[1])
errors1 = convert(Vector{Float64}, random[2])
errors2 = convert(Vector{Float64}, minimized[2])

# input in one vector where first half of vector refers to first group and second half to second group
weights = append!(weights1, weights2)
sx = repeat(["B", "A"], inner = 10)
errors = append!(errors1, errors2)


groupedbar(weights, yerr = errors, group = sx,
        xlabel = L"\textrm{10 bins of width 0.10 from bin 1 = } [0.00,0.10) \textrm{ to bin 10 = } [0.90,1.00)",
        xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, titlefontsize = 10,
        ylabel = L"\textrm{Normalized probability}",
        label = [L"\textrm{Grids with low } G_{av}" L"\textrm{Random grids}"],
        title = L"\textrm{Normalized histogram of flow distribution}",
        bar_width = 0.67,
        lw = 0.0, markerstrokewidth = 0.7, markerstrokecolor = :black,
        c = [:deepskyblue :orange], #https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
        framestyle = :box, grid = true, xticks = 1:1:10)
Plots.savefig("flows_rand_min_Gav.pdf")

## gen-gen, con-con, gen-con
nr_gen_con_av(SData, 3) # random
nr_gen_con_av(SData, 6) # minimized G_av
nr_gen_con_av_diff(SData)

# histogramm of outliers
Data_high_gc_low_Gav = high_gc_low_Gav(SData, 108, 23)
N_side = 10
C = 1.0
Steps_per_temp = 4 # number of steps before lowering temperature
k_max = 10000
ann_shed = "0.999 ^ k"
N_runs = length(Data_high_gc_low_Gav)
T = 0.95
JLD.save("high_gc_low_Gav.jld", "Data", Data_high_gc_low_Gav, "N_side", N_side, "k_max", k_max, "Steps_per_temp", Steps_per_temp, "Threshold T", T, "ann_shed", ann_shed, "N_runs", N_runs, "C", C)
SData_high_gc_low_Gav = load("/media/vollmich/Spass/Runs_Cluster/final_runs3/high_gc_low_Gav.jld")

nbins = 21
random = weights_mean_err(SData, 1, nbins)
minimized = weights_mean_err(SData_high_gc_low_Gav, 4, nbins)

weights1 = convert(Vector{Float64}, random[1])
weights2 = convert(Vector{Float64}, minimized[1])
errors1 = convert(Vector{Float64}, random[2])
errors2 = convert(Vector{Float64}, minimized[2])

# input in one vector where first half of vector refers to first group and second half to second group
weights = append!(weights1, weights2)
sx = repeat(["B", "A"], inner = 20)
errors = append!(errors1, errors2)


groupedbar(weights, yerr = errors, group = sx,
        xlabel = L"\textrm{20 bins of width 0.05 from bin 1: } [0.00,0.05) \textrm{ to bin 20: } [0.95,1.00)",
        xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, titlefontsize = 10,
        ylabel = L"\textrm{Frequency}",
        label = [L"\textrm{Grids with high gc and low } G_{av}" L"\textrm{Random grids}"],
        title = L"\textrm{Distribution: Mean of 1000 random and 17 grids with high gc and low } G_{av}",
        #title = L"\textrm{Flow distribution: Mean of 1000 random and of 1000 grids with minimized } G_{av}",
        bar_width = 0.67,
        lw = 0.0, markerstrokewidth = 0.7, markerstrokecolor = :black,
        c = [:deepskyblue :orange], #https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
        framestyle = :box, grid = true, xticks = 1:1:20)
Plots.savefig("flows_rand_high_gc_low_Gav.pdf")

##locality
locality(SData1)
locality(SData1)

## testing LODF 
g = gen_square_grid(N_side)
g1 = copy(g)
P = gen_stable_config(g, N_side, C)
F1 = flow_old(g1, P)
Flow_first_line = F1[1]
Flow_third_line = F1[3]
visualize_graph(g1, P, N_side)

# our way of recalculating flow after edge removal
rem_edge!(g, 1, 2)
g2 = copy(g)
F2 = flow_old(g2, P)
Flow_third_line = F2[2]
visualize_graph(g2, P, N_side)

Flow_Change = F1 - F2
# LODF-way of recalculating flow after edge removal
l = 1
k = 3

g = gen_square_grid(N_side)
L = Array(LightGraphs.laplacian_matrix(g))
B = Array(incidence_matrix(g, oriented=true))
R = pinv(L)
PTDF = transpose(B) * R
# This is how PTDF is defined: F = PTDF * P

LODF = (PTDF * B)[k][l] / (1 - (PTDF * B)[l][l])
F_new = F1[k] + LODF * F1[l]

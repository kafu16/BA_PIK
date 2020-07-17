using LinearAlgebra

### building of square grid configuration
using LightGraphs, GraphPlot

## generation of non-embedded square grid
function gen_square_grid(N_side) # N_side: this number sqared gives number of vertices, for N_side > 2
    N_vertices = N_side * N_side # N_vertices: number of vertices

    g = LightGraphs.SimpleGraph(N_vertices) # generates undirected graph

    # building of square grid
    # loop one: going through each row of square grid (except last vertex of row, except last row) connecting
    # each vertex with the vertex of its right and with the vertex under it
    for i in 1:N_side:N_vertices - 2 * N_side +1
        for i in i:i + N_side - 2
            LightGraphs.add_edge!(g, i, i+1) # LightGraphs. because of "both Graphs and EmbeddedGraphs export "add_edge!"; uses of it in module Main must be qualified"
            LightGraphs.add_edge!(g, i, i+N_side)
        end
    end

    # loop two: going through last column of square grid connecting each vertex with the vertex under it
    for i in N_side:N_side:N_vertices - N_side
        LightGraphs.add_edge!(g, i, i + N_side)
    end

    # loop three: going through last row of square grid connecting each vertex with the vertex of its right
    for i in N_vertices - N_side + 1:N_vertices -1
        LightGraphs.add_edge!(g, i, i+1)
    end
    g
end


## Building of an array whose entries represent nodes of square grid.
# Vertex one is entry one of array counting until the end of first line
# of square grid, continuing to count in second line and so on.
# Above built graph and array are (state: 25.09.2019) non-interdependent.

#  initial configuration: flow network with random placements nodes having inflow = 1
# or outflow = -1 representing consumers and generators

using Random

function gen_rand_config(N_side) # generates random configuration P (see below)
    N_vertices = N_side * N_side
    P = ones(Float64, N_vertices) # P_i: net inflow or outflow at vertex i
    P[1:2:end] .= P[1:2:end] * -1. # every second element is changed from 1 to -1
    #### ToDo evtl. '.' entfernen, weil nicht einheitlich
    shuffle!(P) # randomly permutes entries of array
end

function gen_stable_config(g, N_side, C) # generation of initial stable random configuration P
#### ToDo build: return error when stable config is not possible due to a too low C
    P = gen_rand_config(N_side)
    F = flow(g, P)
    while maximum(abs.(F)) > C
        P = gen_rand_config(N_side)
        F = flow(g, P)
    end
    P
end


### analysis of square grid configuration

## flow calculation
using IterativeSolvers

# N_vertices: number of vertices
# m: number of edges
function flow(g, P) # calculates the flow of graph g for a certain configuration P
    B = incidence_matrix(g, Float64, oriented=true) # nxm oriented incidence matrix B
    # node-edge matrix, indexed by [v, i], i is in 1:ne(g),
    # indexing edge e. Directed graphs: -1 indicates src(e) == v, 1 indicates dst(e) == v.
    # indexing of the edges e does not correspond to sequence of building edges above!!!
    # indexing of edges: https://github.com/JuliaGraphs/LightGraphs.jl/blob/216f5ffa77860d4c39b8d05fe2197d0af75a4241/src/linalg/spectral.jl#L129-L142

    F = lsqr(B, P) # solving BF = P for F by returning minimal norm solution: https://juliamath.github.io/IterativeSolvers.jl/dev/linear_systems/lsqr/
    F_rounded = round.(F; digits = 2)
end


## failure/removal of edge that causes cascade
function linefailure!(g, i) # edge i is removed
    B = Array(incidence_matrix(g, oriented=true))
    rem_edge!(g, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
    g
end


## cascade
function cascade!(g, P, C) # removes all edges whose flow is bigger than C -> calculates flow again and so on until all flows are smaller than C
    F = flow(g, P)
    while maximum(abs.(F)) > C

        B = Array(incidence_matrix(g, oriented=true))
        m = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
        # it is ne(g) = size(B)[2], where ne() give number of edges
        for i in 1:m # in this loop all edges that carry flow bigger than certain value/ threshold are being removed
            if abs(F[i]) > C
                rem_edge!(g, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
                # B[:, i] gives elements of i-th column of B, isodd accesses only the odd elements of these element
                # findfirst()/ findlast() access first/ last odd element.
                # incidence matrix B comprises only zeros (even elements) and one 1 and one -1 for the two vertices that are connected
                # to respective edge
            end
        end
        F = flow(g, P)
    end
    g
end


## measurements of energy() (see below)
function biggest_component(g) # gives biggest connected component of graph g
    maximum(length.(LightGraphs.connected_components(g))) # connected_components() returns vector whose components contain connected vertices
    # so the length of a component is the size of a connected graph
    # so here the biggest connected graph out of all vertices is chosen
end


## energy function
using Statistics

function energy!(g, P, C, N_side, N_removals = 0) # calculates energy of step k, C: threshold that marks line failure,
# N_side: N_squared gives number of vertices, N_removals (optional argument): number of edge removals for approximation
    B = Array(incidence_matrix(g, oriented=true))
    m = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
    linefailure_indizes = collect(1:m) # collect() collects all values in the range 1:m in an array, here all edges are numberd in an array
    linefailure_indizes = shuffle!(linefailure_indizes) # positions of edges in linefailure_indizes is randomly permuted by shuffle!()
    # randomness is redundant unsless m in following for-loop is replaced by a number smaller than m

    if N_removals > 0
        N = N_removals
    else
        N = m
    end

    G = zeros(N)

    #### ToDo Für Näherung, d.h. man zieht nur N edges anstatt alle, muss man das nachfolgende m durch N ersetzen -> herausfinden ob man das iwie mit optional arguments machen kann
    for i in 1:N # for loop for randomly chosen N_vertices linefailures
        g = linefailure!(g, linefailure_indizes[i])
        g = cascade!(g, P, C)
        #global G[i] = ne(g)
        global G[i] = biggest_component(g) # G: size of biggest connected component, global lets the values of G saved after each run of loop,
        # otherwise G would be overwritten in each run of loop
        g = gen_square_grid(N_side)
    end
    G_av = mean(G)
    G, G_av # this way two values in a tuple are returned by a function
end

## swap of configurations
function swap!(P, N_side) # swaps randomly chosen generator-consumer pair
    N_vertices = N_side * N_side
    A = rand(1:N_vertices,1)
    B = rand(1:N_vertices,1)
    # println("A initial", A) # following println-lines only to check what while-loop does, not necessary for swap()
    # println("B initial", B)
    while A == B || P[A] == P[B]
    # A == B avoids that entity is exchanged by itself, P[A] = P[B] avoids consumer-consumer and generator-generator-swap,
    # while-loop is run if condition on either left or right hand side of || is fulfilled
        # println("A while condition", A, "P[A] while condition", P[A])
        # println("B while condition", B, "P[B] while condition", P[B])
        A = rand(1:N_vertices,1) # if one writes "!" after function, behaviour of function changes see https://docs.julialang.org/en/v1/stdlib/Random/
        B = rand(1:N_vertices,1) #### ToDo warum Fehlermeldung, wenn man A und B als global assigned?
        # println("A in while", A)
        # println("B in while", B)
    end
    # println("A after while", A, "P[A] after while", P[A])
    # println("B after while", B, "P[B] after while", P[B])
    P[A], P[B] = P[B], P[A] # swaps neighbors
    # println(P)
    P
end


# generation of stable swapped configuration
function stable_swapped_config!(g, P, C, N_side)
# to avoid calculation steps, the input configuration should be stable as the stable_swapped_config() only permutes
# one generator-consumer pair. so having an unstable configuration as input will probably take more steps than
# first generating a stable configuration by gen_stable_config() and then apply stable_swapped_config()
#### ToDo build: return error when stable config is not possible due to a too low C
    P_stable_old = copy(P)
    # by using P_stable_old it is made sure that a the swapped configuration differs only by one permutation
    # otherwise by the following while-loop in each iteration of the loop the newly generated configuration would be swapped
    # again to obtain a stable configuration. This can not be done by P_stable_old = P because of variable assignment. In the latter case
    # P_stable_old would equal P all the time, even when the value of P is changed
    P = swap!(P, N_side)
    F = flow(g, P)
    while maximum(abs.(F)) > C
        P = swap!(P_stable_old, N_side)
        F = flow(g, P)
    end
    P
end


## evaluation of energy one random configuration
N_removals = 0
# """
# function for evaluating the biggest component without evaluating several other functions
# """
function eval_one_random_config(N_side, C)
    g = gen_square_grid(N_side)
    P = gen_stable_config(g, N_side, C)
    energy!(g, P, C, N_side, N_removals)
end



### Embedding of Simulated Annealing (SA)

## key code of SA

#### ToDo: Option N_removals entfernen

# core function for simulated annealing
function sim_anneal!(g, P, C, N_side, N_removals = 0, k_max = 10) # k_max: setting number of computation steps
    # given an initial configuration P sim_anneal!() tendentially finds a more stable configuration
    en = [ ]
    for k in 0:k_max - 1
        T = temperature(k) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
        P_old = copy(P) # for calculating the energy of "old" configuration
        P = stable_swapped_config!(g, P, C, N_side)
        energy_old = energy!(g, P_old, C, N_side, N_removals)[2] # by [2] only the second value of tuple is returned (G_av)
        g = gen_square_grid(N_side) # energy!() mutates g, so g has to be rebuilt every time before calculating energy!() again
        energy_new = energy!(g, P, C, N_side, N_removals)[2] # by [2] only the second value of tuple is returned (G_av)
        ΔE = energy_new - energy_old
        #### performance: let energy() calculate G_av only

        if ΔE <= 0 # man könnte conditional auch umdrehen: if (ΔE <= 0 AND probability(ΔE, T) < rand())
                                                                 # P = P_old
            P
        elseif probability(ΔE, T) > rand() # rand() gives random number element of [0,1]
            P
        else
            P = P_old
        end
        g = gen_square_grid(N_side)
        push!(en, energy_old)
    end
    P, en
end

# applies SA on random square grid
#### ToDo evtl. diesen code direkt in sim_anneal() einbauen. Bin mir unsicher, ob das sinnvoll ist...
function eval_sim_anneal!(N_side, C, T, N_removals = 0, k_max = 10)
    g = gen_square_grid(N_side)
    P = gen_stable_config(g, N_side, C) # to avoid iteration steps it is important to start with a stable configurations see comment at stable_swapped_config!()
    P_initial = copy(P)
    P, en = sim_anneal!(g, P, C, N_side, 0, k_max)
    g = gen_square_grid(N_side)
    energy_initial = energy!(g, P_initial, C, N_side, N_removals)
    g = gen_square_grid(N_side)
    N_T = flows_above_thres(T, P, g)
    energy_final = energy!(g, P, C, N_side, N_removals)
    P_initial, energy_initial, P, energy_final, N_T, en
end

function collect_data_SA_runs_var_ann_shed(N_runs, N_side, C, T, annealing_schedule, k_max)
    Data = []
    for i in 1:N_runs
        g = gen_square_grid(N_side)
        P = gen_stable_config(g, N_side, C) # to avoid iteration steps it is important to start with a stable configurations see comment at stable_swapped_config!()
        P_initial = copy(P)
        en = [ ]
        N_removals = 0
        for k in 0:k_max - 1
            Temp = annealing_schedule(k)
            P_old = copy(P) # for calculating the energy of "old" configuration
            P = stable_swapped_config!(g, P, C, N_side)
            energy_old = energy!(g, P_old, C, N_side, N_removals)[2] # by [2] only the second value of tuple is returned (G_av)
            g = gen_square_grid(N_side) # energy!() mutates g, so g has to be rebuilt every time before calculating energy!() again
            energy_new = energy!(g, P, C, N_side, N_removals)[2] # by [2] only the second value of tuple is returned (G_av)
            ΔE = energy_new - energy_old
            #### performance: let energy() calculate G_av only

            if ΔE <= 0 # man könnte conditional auch umdrehen: if (ΔE <= 0 AND probability(ΔE, T) < rand())
                                                                     # P = P_old
                P
            elseif probability(ΔE, Temp) > rand() # rand() gives random number element of [0,1]
                P
            else
                P = P_old
            end
            g = gen_square_grid(N_side)
            push!(en, energy_old)
        end
        g = gen_square_grid(N_side)
        energy_initial = energy!(g, P_initial, C, N_side, N_removals)
        g = gen_square_grid(N_side)
        N_T = flows_above_thres(T, P_initial, g), flows_above_thres(T, P, g)
        locality_init = loc_1step!(P_initial, C, N_side), loc_1step_0!(P_initial, C, N_side)
        g = gen_square_grid(N_side)
        locality_final = loc_1step!(P, C, N_side), loc_1step_0!(P, C, N_side)
        g = gen_square_grid(N_side)
        energy_final = energy!(g, P, C, N_side, N_removals)
        SA_extremal = P_initial, energy_initial, nr_gen_con(P_initial, N_side), P, energy_final, nr_gen_con(P, N_side), N_T, en, locality_init, locality_final
        push!(Data, SA_extremal)
    end
    Data
end


## several functions for SA

## annealing schedules
function temperature(k) # arbitrary temperature function that decreases to zero and calculates temperature dependent of step
    1. / (1 + 0.0001 * floor(k/4))
end

function temp1(k) # arbitrary temperature function that decreases to zero and calculates temperature dependent of step
    1. / (1 + 0.0001 * floor(k/4))
end

function temp2(k)
    1. / (1 + 0.0007 * floor(k/4))
end

function temp3(k)
    1. / (1 + 0.0015 * floor(k/4))
end

function temp4(k)
    1. / (1 + 0.0025 * floor(k/4))
end

function temp5(k)
    10. / (1 + 0.02 * k)
end

function temp_ex1(k)
    0.999 ^ (floor(k/4)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex2(k)
    0.9995 ^ (floor(k/4)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex3(k)
    0.99975 ^ (floor(k/4)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex4(k)
    0.9999 ^ (floor(k/4)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end



function probability(ΔE, T) # probability function depends on ΔE and on temperature function
    exp( - ΔE / T) # Boltzmann's constant k_B is set to one
end


### visualisation of square grid

## embed square grid
function set_vertex_locs(N_side) # sets vertex locations, example usage: gplot(g, locs_x, locs_y) while locs_x, locs_y = set_vertex_locs()
    N_vertices = N_side * N_side
    locs_x = zeros(Int64, N_vertices) # generation of arrays
    locs_y = zeros(Int64, N_vertices)

    # generation of locs_x
    k = 1
    j = 1
    for i in 1:N_side
        for j in 1:N_side
            locs_x[k] = j
            k = k + 1
            j = j + 1
        end
        j = 1
    end

    # generation of locs_y
    k = 1
    for i in 1:N_side
        for j in 1:N_side
            locs_y[k] = i
            k = k + 1
        end
    end

    locs_x = convert(Vector{Float64}, locs_x) # conversion from Vector{Int64} to Vector{Float64}
    locs_y = convert(Vector{Float64}, locs_y)
    locs_x, locs_y
end

## show P_i in Graph
using Graphs

# P = gen_rand_config(N_side)
# nodelabel = P
# gplot(g, nodelabel = P)

## colour P_i
using Colors

# choice of colours: https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
function set_colours(P)
    nodefillc = [] # generates empty vector of type 'Any'
    for value in P
        if value == 1
            push!(nodefillc, colorant"grey85") # push! inserts items at end of collection # in TeX grey 95
        else
            push!(nodefillc, colorant"black") # in TeX grey 75
        end
    end
    nodefillc
end

function set_colours2(P)
    nodefillc = [] # generates empty vector of type 'Any'
    for value in P
        if value == 1
            push!(nodefillc, colorant"grey95") # push! inserts items at end of collection # in TeX grey 95
        else
            push!(nodefillc, colorant"grey75") # in TeX grey 75
        end
    end
    nodefillc
end
## show Flows F_e in Graph
# F = flow(g, P)
# g = gplot(g, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)

## visualization of graph
function visualize_graph(g, P, N_side)
    F = flow(g, P)
    locs_x, locs_y = set_vertex_locs(N_side)
    nodefillc = set_colours(P)

    gplot(g, locs_x, locs_y, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
    #gplot(g, locs_x, locs_y, nodelabel = P, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
    #gplot(g, locs_x, locs_y, nodefillc = nodefillc)
end

function visualize_graph_vlabel(g, P, N_side)
    F = flow(g, P)
    locs_x, locs_y = set_vertex_locs(N_side)
    nodefillc = set_colours2(P)

    #gplot(g, locs_x, locs_y, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
    gplot(g, locs_x, locs_y, nodelabel = P, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
    #gplot(g, locs_x, locs_y, nodefillc = nodefillc)
end

## flow direction, definition of P_i:
## (plotted SimpleDiGraph: all horizontal lines point to the right, all vertical line point downwards)
# if F_e < 0, flow rightwards or downwards (flow is in the same direction as arrow)
# if F_e > 0, flow leftwards or upwards (flow is in the opposite direction as arrow)
# P_i = 1 one unit of flow is generated
# P_i = -1 one unit of flow is consumed

## visualize graph after line failure induced cascade
# for evaluation of visualize_graph_after_linefailure_cascade: it must be line ⋜ m:
#B = Array(incidence_matrix(g, oriented=true))
#m = size(B)[2]
function visualize_graph_after_linefailure_cascade(P, C, N_side, line)
    g = gen_square_grid(N_side)
    g = linefailure!(g, line)
    g = cascade!(g, P, C)
    visualize_graph(g, P, N_side)
end

### data colleciton

## data collection: collects multiple runs of eval_sim_anneal!() in one object
function collect_data_SA_runs(N_runs, N_side, C, T, k_max)
    Data = []
    for i in 1:N_runs
        SA_extremal = eval_sim_anneal!(N_side, C, T, 0, k_max)
        push!(Data, SA_extremal)
    end
    Data
end

## safe data
using JLD
using Dates #### ToDo for saving date and time in filename
#Dates.now(Dates.UTC)

## visualize data
function visualize_data(SData, rand_opt, Run_Nr) # for random grid: rand_opt=1 for optimized grid: rand_opt=4
    Data = SData["Data"]
    N_side = SData["N_side"]
    P = Data[Run_Nr][rand_opt]
    g = gen_square_grid(N_side)
    visualize_graph(g, P, N_side)
end

function flows_above_thres(T, P_SA, g) # gives number of flows above threshold T
    F = flow(g, P_SA)
    count(x -> x > T, abs.(F))
end

function collect_data(Data, col)
    Data_x = [ ]
    N_runs = length(Data)
    for i in 1:N_runs
        x = Data[i][col]
        push!(Data_x, x)
    end
    Data_x
end

function collect_data2(Data, col, subcol)
    Data_x = [ ]
    N_runs = length(Data)
    for i in 1:N_runs
        x = Data[i][col][subcol]
        push!(Data_x, x)
    end
    Data_x
end

function collect_av_e_std(SData)
    Data = SData["Data"]
    N_runs = SData["N_runs"]
    Data_mean = [ ]
    Data_std = [ ]
    N_steps = length(Data[1][8])
    for i in 1:N_steps
        x = mean(collect_data2(Data, 8, i))
        push!(Data_mean, x)
        y = 1 / sqrt(N_runs) * std(collect_data2(Data, 8, i); corrected=true)
        push!(Data_std, y)
    end
    Data_mean, Data_std
end

function collect_av_e_diff_std(Data)
    Data_mean = [ ]
    Data_std = [ ]
    N_runs = length(Data[1][6])
    for i in 2:N_runs
        x = mean(collect_data2(Data, 6, i - 1)) - mean(collect_data2(Data, 6, i))
        push!(Data_mean, x)
        y = sqrt((1 / sqrt(N_runs) * std(collect_data2(Data, 6, i - 1); corrected=true))^2 + (1 / sqrt(N_runs) * std(collect_data2(Data, 6, i); corrected=true))^2)
        push!(Data_std, y)
    end
    Data_mean, Data_std
end

function energy_from_data(Data, N_runs, C, N_side) # davor energy() abändern
    edge_energy = [ ]
    for i in 1:N_runs
        g = gen_square_grid(N_side)
        energy = energy!(g, Data[i][3], C, N_side)
        push!(edge_energy, energy[2])
    end
    #append!(Data, edge_energy)
    #Data
    edge_energy
end

function energy_from_data2(Data, N_runs, C, N_side)
    g_av_energy = [ ]
    for i in 1:N_runs
        energy = Data[i][4][2]
        push!(g_av_energy, energy)
    end
    #append!(Data, edge_energy)
    #Data
    g_av_energy
end


# histograms
function flow_single_sample(SData, i, j) # i: sample number, j = 1: random, j = 4: optimized
    N_side = SData["N_side"]
    g = gen_square_grid(N_side)
    Data = SData["Data"]
    P = Data[i][j]
    F = flow(g, P)
    x = abs.(F)
end

function flows_N_runs(SData, j) # i: sample number, j = 1: random grids, j = 4: optimised grids
    N_side = SData["N_side"]
    g = gen_square_grid(N_side)
    Data = SData["Data"]
    All_flows = [ ]
    N_runs = length(Data)
    P = Data[1][j]
    F = flow(g, P)
    All_flows = copy(F)
    for i in 2:N_runs
        P = Data[i][j]
        F = flow(g, P)
        append!(All_flows, F)
    end
    All_flows
    x = abs.(All_flows)
end

function high_gc_low_Gav(SData, gen_con, G_av_final)
    Data = SData["Data"]
    N_runs = length(Data)
    configs_high_gc_low_Gav = [ ]
    for i in 1:N_runs
        if Data[i][6][3] > gen_con && Data[i][5][2] < G_av_final # 23 is mean G_av plus STD
            configs_high_gc_low_Gav = push!(configs_high_gc_low_Gav, Data[i])
        end
    end
    configs_high_gc_low_Gav
end

function flows_high_gc_low_Gav(Data, N_side, j) # i: sample number, j = 1: random grids, j = 4: optimised grids
    g = gen_square_grid(N_side)
    All_flows = [ ]
    N_runs = length(Data)
    P = Data[1][j]
    F = flow(g, P)
    All_flows = copy(F)
    for i in 2:N_runs
        P = Data[i][j]
        F = flow(g, P)
        append!(All_flows, F)
    end
    All_flows
    x = abs.(All_flows)
end

function weights_mean_err(SData, j, nbins) # j = 1: random grids, j = 4: optimised grids
    weights = [ ]
    N_bars = nbins - 1
    for i in 1:N_bars
        push!(weights, [ ])
    end
    Data = SData["Data"]
    N_runs = length(Data)
    bins = range(0.0, 1.0, length = nbins)
    for i in 1:N_runs
        flows = flow_single_sample(SData, i, j)
        h = fit(Histogram, flows, bins)
        h = normalize(h, mode=:probability)
        weightvalues = h.weights
        for i in 1:N_bars
            append!(weights[i], weightvalues[i])
        end
    end
    weights_av = [ ]
    weights_err = [ ]
    for i in 1:N_bars
        append!(weights_av, mean(weights[i]))
        append!(weights_err, 1 / sqrt(N_runs) * Statistics.std(weights[i]))
    end
    weights_av, weights_err
end

# measure  number of edges gen-gen, con-con, gen-con 04.05.2020
function nr_gen_con(P, N_side)
    N_vertices = N_side * N_side
    gen_gen = 0 # generates empty vector of type 'Any'
    con_con = 0
    gen_con = 0
    for i in 1:N_side:N_vertices - N_side + 1
        for i in i:i + N_side - 2
            if P[i] == 1 && P[i + 1] == 1
                gen_gen = gen_gen + 1
            elseif P[i] == -1 && P[i + 1] == -1
                con_con = con_con + 1
            else gen_con = gen_con + 1
            end
        end
    end
    for i in 1:N_side
        for i in i:N_side:N_vertices - (2 * N_side) + i
            if P[i] == 1 && P[i + N_side] == 1
                gen_gen = gen_gen + 1
            elseif P[i] == -1 && P[i + N_side] == -1
                con_con = con_con + 1
            else gen_con = gen_con + 1
            end
        end
    end
    gen_gen, con_con, gen_con
end

function nr_gen_con_av(SData, col)
    Data = SData["Data"]
    N_runs = SData["N_runs"]
    gen_gen_av = mean(collect_data2(Data, col, 1))
    con_con_av = mean(collect_data2(Data, col, 2))
    gen_con_av = mean(collect_data2(Data, col, 3))
    gen_gen_std = 1 / sqrt(N_runs) * std(collect_data2(Data, col, 1))
    con_con_std = 1 / sqrt(N_runs) * std(collect_data2(Data, col, 2))
    gen_con_std = 1 / sqrt(N_runs) * std(collect_data2(Data, col, 3))
    gen_gen_av, gen_gen_std, con_con_av, con_con_std, gen_con_av, gen_con_std
end

function nr_gen_con_av_diff(SData)
    Data = SData["Data"]
    N_runs = SData["N_runs"]
    gen_gen_av = mean(collect_data2(Data, 3, 1) - collect_data2(Data, 6, 1))
    con_con_av = mean(collect_data2(Data, 3, 2) - collect_data2(Data, 6, 2))
    gen_con_av = mean(collect_data2(Data, 3, 3) - collect_data2(Data, 6, 3))
    gen_gen_std = sqrt((1 / sqrt(N_runs) * std(collect_data2(Data, 3, 1)))^2 + (1 / sqrt(N_runs) * std(collect_data2(Data, 6, 1)))^2) # the measurement uncertainties sum up
    con_con_std = sqrt((1 / sqrt(N_runs) * std(collect_data2(Data, 3, 2)))^2 + (1 / sqrt(N_runs) * std(collect_data2(Data, 6, 2)))^2)
    gen_con_std = sqrt((1 / sqrt(N_runs) * std(collect_data2(Data, 3, 3)))^2 + (1 / sqrt(N_runs) * std(collect_data2(Data, 6, 3)))^2)
    gen_gen_av, gen_gen_std, con_con_av, con_con_std, gen_con_av, gen_con_std
end


# locality functions
function cascading_steps!(g, P, C) # removes all edges whose flow is bigger than C -> calculates flow again and so on until all flows are smaller than C
    F = flow(g, P)
    cascading_steps = [ ]
    while maximum(abs.(F)) > C

        B = Array(incidence_matrix(g, oriented=true))
        m = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
        # it is ne(g) = size(B)[2], where ne() give number of edges
        for i in 1:m # in this loop all edges that carry flow bigger than certain value/ threshold are being removed
            if abs(F[i]) > C
                rem_edge!(g, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
                # B[:, i] gives elements of i-th column of B, isodd accesses only the odd elements of these element
                # findfirst()/ findlast() access first/ last odd element.
                # incidence matrix B comprises only zeros (even elements) and one 1 and one -1 for the two vertices that are connected
                # to respective edge
            end
        end
        g_new = copy(g)
        push!(cascading_steps, g_new)
        F = flow(g, P)
    end
    cascading_steps
end

function cascade_1step!(g, P, C) # removes all edges whose flow is bigger than C, outputs vertices of those edges
    F = flow(g, P)
    if maximum(abs.(F)) > C
        B = Array(incidence_matrix(g, oriented=true))
        n = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
        # it is ne(g) = size(B)[2], where ne() give number of edges
        vertices_of_failed_edges = [ ]
        for i in 1:n # in this loop all edges that carry flow bigger than certain value/ threshold are being removed
            if abs(F[i]) > C
                rem_edge!(g, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
                # B[:, i] gives elements of i-th column of B, isodd accesses only the odd elements of these element
                # findfirst()/ findlast() access first/ last odd element.
                # incidence matrix B comprises only zeros (even elements) and one 1 and one -1 for the two vertices that are connected
                # to respective edge
                push!(vertices_of_failed_edges, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
            end
        end
    else vertices_of_failed_edges = [ ]
    end
    vertices_of_failed_edges
end

function loc_1step!(P, C, N_side) # calculates average distance of linefailure from initial line failure for all edges
# N_side: N_squared gives number of vertices, N_removals (optional argument): number of edge removals for approximation
    g = gen_square_grid(N_side)
    B = Array(incidence_matrix(g, oriented=true))
    n = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
    linefailure_indizes = collect(1:n) # collect() collects all values in the range 1:m in an array, here all edges are numberd in an array
    #linefailure_indizes = shuffle!(linefailure_indizes) # positions of edges in linefailure_indizes is randomly permuted by shuffle!()

    min_failure_distances = [ ]
    for i in 1:n
        g = gen_square_grid(N_side)
        B = Array(incidence_matrix(g, oriented=true))
        first = findfirst(isodd, B[:, i])
        last = findlast(isodd, B[:, i])
        rem_edge!(g, first, last)
        vertices_initial_failed_edge = [first, last]
        vertices_of_failed_edges = cascade_1step!(g, P, C)
        if isempty(vertices_of_failed_edges)
            continue
        end
        failure_distances = [ ]
        for i in 1:length(vertices_initial_failed_edge)
            for j in 1:length(vertices_of_failed_edges)
                push!(failure_distances, length(a_star(g,vertices_initial_failed_edge[i],vertices_of_failed_edges[j])))
            end
        end
        push!(min_failure_distances, minimum(failure_distances))
    end
    if isempty(min_failure_distances)
        min_failure_distances = [0]
    end
    min_failure_distance_av = mean(min_failure_distances)
    min_failure_distances, min_failure_distance_av
end

function loc_1step_0!(P, C, N_side) # only differing to loc_1step!() by counting the distance
# as zero in case no other edge failure is caused by initial edge removal
    g = gen_square_grid(N_side)
    B = Array(incidence_matrix(g, oriented=true))
    n = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
    linefailure_indizes = collect(1:n) # collect() collects all values in the range 1:m in an array, here all edges are numberd in an array
    #linefailure_indizes = shuffle!(linefailure_indizes) # positions of edges in linefailure_indizes is randomly permuted by shuffle!()

    min_failure_distances = [ ]
    for i in 1:n
        g = gen_square_grid(N_side)
        B = Array(incidence_matrix(g, oriented=true))
        first = findfirst(isodd, B[:, i])
        last = findlast(isodd, B[:, i])
        rem_edge!(g, first, last)
        vertices_initial_failed_edge = [first, last]
        vertices_of_failed_edges = cascade_1step!(g, P, C)
        failure_distances = [ ]
        if isempty(vertices_of_failed_edges)
            push!(failure_distances, 0)
        end

        for i in 1:length(vertices_initial_failed_edge)
            for j in 1:length(vertices_of_failed_edges)
                push!(failure_distances, length(a_star(g,vertices_initial_failed_edge[i],vertices_of_failed_edges[j])))
            end
        end
        push!(min_failure_distances, minimum(failure_distances))
    end
    min_failure_distance_av = mean(min_failure_distances)
    min_failure_distances, min_failure_distance_av
end

function locality(SData)
    Data = SData["Data"]
    N_side = SData["N_side"]
    C = SData["C"]
    N_runs = SData["N_runs"]
    configs_init= collect_data(Data, 1)
    configs_final = collect_data(Data, 4)
    Data_init = [ ]
    Data_final = [ ]
    Data_init0 = [ ]
    Data_final0 = [ ]
    N_runs = length(Data)
    for i in 1:N_runs
        x = loc_1step!(configs_init[i], C, N_side)
        push!(Data_init, x)
        y = loc_1step!(configs_final[i], C, N_side)
        push!(Data_final, y)
        v = loc_1step_0!(configs_init[i], C, N_side)
        push!(Data_init0, v)
        w = loc_1step_0!(configs_final[i], C, N_side)
        push!(Data_final0, w)
    end
    locality_av = mean(collect_data(Data_final, 2) - collect_data(Data_init, 2))
    locality_std = sqrt((1 / sqrt(N_runs) * std(collect_data(Data_final, 2)))^2 + (1 / sqrt(N_runs) * std(collect_data(Data_init, 2)))^2)
    locality_av0 = mean(collect_data(Data_final0, 2) - collect_data(Data_init0, 2))
    locality_std0 = sqrt((1 / sqrt(N_runs) * std(collect_data(Data_final0, 2)))^2 + (1 / sqrt(N_runs) * std(collect_data(Data_init0, 2)))^2)
    #locality_std0 = 1 / sqrt(N_runs) * (std(collect_data(Data_final0, 2)) + std(collect_data(Data_init0, 2)))
    A = locality_av, locality_std, Data_init, Data_final, mean(collect_data(Data_init, 2)), 1 / sqrt(N_runs) * std(collect_data(Data_init, 2)), mean(collect_data(Data_final, 2)), 1 / sqrt(N_runs) * std(collect_data(Data_final, 2))
    B = locality_av0, locality_std0, Data_init0, Data_final0, mean(collect_data(Data_init0, 2)), 1 / sqrt(N_runs) * std(collect_data(Data_final0, 2)), mean(collect_data(Data_final0, 2)), 1 / sqrt(N_runs) * std(collect_data(Data_init0, 2))
    A, B
end


using Plots
## export of plots
using Compose
using Cairo
using Fontconfig
#
# # save to pdf
# draw(PDF("bla.pdf", 16cm, 16cm), gplot(g))
# # save to png
# draw(PNG("bla.png", 16cm, 16cm), gplot(g))
# # save to svg
# draw(SVG("bla_neu.svg", 16cm, 16cm), gplot(g))

# performance
using CPUTime

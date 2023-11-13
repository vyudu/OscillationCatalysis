using Graphs, MetaGraphs
using LinearAlgebra, Combinatorics

function butlerVolmer(E, E_eq, C_O, C_R; α=0.5, β=0.0004036, j_0=1.)
    f = 96500 * β
    return j_0*exp(-α*f*(E - E_eq)), j_0*exp((1-α)*f*(E - E_eq))
end

function matrixTree(g::MetaDiGraph)
    # generate all spanning trees
    n = nv(g)
    ug = SimpleGraph(SimpleDiGraph(g))
    trees = collect(combinations(collect(edges(ug)), n-1))
    trees = SimpleGraph.(trees)
    trees = filter!(t->isempty(cycle_basis(t)), trees)
    
    π_i = Float64[]
    W(t) = prod([get_prop(g, e, :rate) for e in edges(t)])

    # constructed rooted trees for every edge
    for v in 1:n
        rootedTrees = [reverse(bfs_tree(t, v, dir=:in)) for t in trees]
        push!(π_i, sum([W(t) for t in rootedTrees]))
    end
    
    # sum the contributions
    return π_i ./ sum(π_i)
end

function steadyState(R::Matrix{Float64}) 
    v = eigvecs(R)
    real.(v[:,end] ./ sum(v[:,end]))
end

### Parameters
# E = electrode potential
# E_eq = equilibrium potential of reaction
# n_e = number of electrons transferred
# C = concentrations of O, R
# R0 = intrinsic rates of downstream reactions
# ΔG = reaction barriers of downstream reactions
# β = inverse temperature
# α = asymmetry factors

function rateMatrix(E, E_eq::Vector{Float64}, n_e::Int64, C::Vector{Float64}; R0 = [1., 10.], ΔG = [1.,2.], β = 0.0004036, α=[0.5,0.5]) 
    n_states = n_e + 3
    R = zeros(n_states, n_states)
    # set rates for Butler-Volmer transitions
    for i in 1:n_e
        R[i+1,i], R[i,i+1] = butlerVolmer(E, E_eq[i], C..., β = β, α=α[i]) 
    end
    # set rates for non-BV transitions - loop 1
    R[n_e+1,n_e+2], R[n_e+2,n_e+1] = R0[1], R0[1]*exp(-ΔG[1])
    R[n_e+2,1], R[1,n_e+2] = R0[1], R0[1]*exp(-ΔG[1])
    # Loop 2
    R[n_e+1,n_e+3], R[n_e+3,n_e+1] = R0[2], R0[2]*exp(-ΔG[2])
    R[n_e+3,1], R[1,n_e+3] = R0[2], R0[2]*exp(-ΔG[2])

    for i in 1:n_states
        totalRate = sum(R[:,i])
        R[i,i] = -totalRate
    end
    
    # Construct graph
    graph = MetaDiGraph(SimpleDiGraph(R))
    for i in CartesianIndices(R)
        set_prop!(graph, i[2], i[1], :rate, R[i])
    end
    return R, graph
end

function selectivityPlot(E_l, E_h, E_eq::Vector{Float64}; α = [0.5, 0.5])
    E_range = LinRange(E_l, E_h, 100)
    selectivityArr = Float64[]
    for E in E_range
        R, g = rateMatrix(E, E_eq, 2, [1., 1.], α = α)
        ss = matrixTree(g)
        s = selectivity(R, ss)
        push!(selectivityArr, s)
    end
    p = plot(E_range, selectivityArr, title="Selectivity vs. Electrode Voltage", xlabel="Voltage", ylabel="Selectivity")
    p
end

function selectivity(R, ss; i=3)
    (R[i+1,i]*ss[i] - R[i,i+1]*ss[i+1]) / (R[i+2,i]*ss[i] - R[i,i+2]*ss[i+2])
end

function affinity(R)
    log(R[2,1] * R[3,2] * R[4,3] * R[1,4] / (R[1,2] * R[2,3] * R[3,4] * R[4,1]))
end

function current1(R, ss; i=3)
    R[i+1,i]*ss[i] - R[i,i+1]*ss[i+1]
end

function current2(R, ss; i=3)
    R[i+2,i]*ss[i] - R[i,i+2]*ss[i+2]
end

function oscillationAffinity(ΔE, E_eq::Vector{Float64}; df = 0.5, E0 = 0.0, α = [0.5, 0.5])
    ΔE_range = LinRange(E0, E0 + ΔE, 100)
    affinity_l = Float64[]
    affinity_h = Float64[]
    affinity_a = Float64[]
    
    for dE in ΔE_range
        E_a = (E_eq[1] + E_eq[2]) / 2 ; E_l = E_a - dE; E_h = E_a + dE
        R_l, g_l = rateMatrix(E_l, E_eq, 2, [1., 1.], α = α)
        R_h, g_h = rateMatrix(E_h, E_eq, 2, [1., 1.], α = α)
        R_a = df*R_l + (1-df)*R_h
        a_l, a_h, a_a = affinity(R_l), affinity(R_h), affinity(R_a)        
        push!(affinity_l, a_l); push!(affinity_h, a_h); push!(affinity_a, a_a)
    end
    p = plot(ΔE_range, [affinity_l, affinity_h, affinity_a], label=["Low" "High" "Oscillating"], title="Affinity vs. Oscillation Magnitude", xlabel="Oscillation Magnitude from Average E", ylabel="Affinity")
    p
end

function oscillationCurrent(ΔE, E_eq::Vector{Float64}; df = 0.5, E0 = 0.0, α = [0.5, 0.5])
    ΔE_range = LinRange(E0, E0 + ΔE, 100)
    current_l1 = Float64[]
    current_h1 = Float64[]
    current_a1 = Float64[]
    current_l2 = Float64[]
    current_h2 = Float64[]
    current_a2 = Float64[]
    
    for dE in ΔE_range
        E_a = (E_eq[1] + E_eq[2]) / 2 ; E_l = E_a - dE; E_h = E_a + dE
        R_l, g_l = rateMatrix(E_l, E_eq, 2, [1., 1.], α = α)
        R_h, g_h = rateMatrix(E_h, E_eq, 2, [1., 1.], α = α)
        R_a = df*R_l + (1-df)*R_h
        ss_l, ss_h, ss_a = matrixTree(g_l), matrixTree(g_h), steadyState(R_a)
        c_l1, c_h1, c_a1 = current1(R_l, ss_l), current1(R_h, ss_h), current1(R_a, ss_a)        
        push!(current_l1, c_l1); push!(current_h1, c_h1); push!(current_a1, c_a1)
        
        c_l2, c_h2, c_a2 = current2(R_l, ss_l), current2(R_h, ss_h), current2(R_a, ss_a)        
        push!(current_l2, c_l2); push!(current_h2, c_h2); push!(current_a2, c_a2)
    end
    p = plot(ΔE_range, [current_l1, current_h1, current_a1, current_l2, current_h2, current_a2], label=["Loop 1 (Low)" "Loop 1 (High)" "Loop 1 (Oscillating)" "Loop 2 (Low)" "Loop 2 (High)" "Loop 2 (Oscillating)"], title="Current vs. Oscillation Magnitude", xlabel="Oscillation Magnitude from Average E", ylabel="Current")
    p
end

function oscillationSelectivity(ΔE, E_eq::Vector{Float64}; df = 0.5, E0 = 0.0, α = [0.5, 0.5])
    ΔE_range = LinRange(E0, E0 + ΔE, 100)
    selectivity_l = Float64[]
    selectivity_h = Float64[]
    selectivity_a = Float64[]
    
    for dE in ΔE_range
        E_a = (E_eq[1] + E_eq[2]) / 2 ; E_l = E_a - dE; E_h = E_a + dE
        R_l, g_l = rateMatrix(E_l, E_eq, 2, [1., 1.], α = α)
        R_h, g_h = rateMatrix(E_h, E_eq, 2, [1., 1.], α = α)
        R_a = df*R_l + (1-df)*R_h
        ss_l, ss_h, ss_a = matrixTree(g_l), matrixTree(g_h), steadyState(R_a)
        s_l, s_h, s_a = selectivity(R_l, ss_l), selectivity(R_h, ss_h), selectivity(R_a, ss_a)
        push!(selectivity_l, s_l); push!(selectivity_h, s_h); push!(selectivity_a, s_a)
    end
    p = plot(ΔE_range, [selectivity_l, selectivity_h, selectivity_a], label=["Low" "High" "Oscillating"], title="Selectivity vs. Oscillation Magnitude", xlabel="Oscillation Magnitude from Average E", ylabel="Selectivity")
    p
end

#TODO:
function oscillationSelectivity3D(ΔE, E_eq::Vector{Float64}; α = [0.5, 0.5])
    ΔE_range = LinRange(0, ΔE, 100)
    selectivity_l = Float64[]
    selectivity_h = Float64[]
    selectivity_a = Float64[]
    
    for dE in ΔE_range
        E_a = (E_eq[1] + E_eq[2]) / 2 ; E_l = E_a - dE; E_h = E_a + dE
        R_l, g_l = rateMatrix(E_l, E_eq, 2, [1., 1.], α = α)
        R_h, g_h = rateMatrix(E_h, E_eq, 2, [1., 1.], α = α)
        R_a = (R_l + R_h) / 2
        ss_l, ss_h, ss_a = matrixTree(g_l), matrixTree(g_h), steadyState(R_a)
        s_l, s_h, s_a = selectivity(R_l, ss_l), selectivity(R_h, ss_h), selectivity(R_a, ss_a)
        push!(selectivity_l, s_l); push!(selectivity_h, s_h); push!(selectivity_a, s_a)
    end
    p = plot(ΔE_range, [selectivity_l, selectivity_h, selectivity_a], label=["Low" "High" "Oscillating"], title="Selectivity vs. Oscillation Magnitude", xlabel="Oscillation Magnitude from Average E", ylabel="Selectivity")
    p
end

# TODO: Oscillation selectivity plots
# Current plots 


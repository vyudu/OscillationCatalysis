#Tunable parameters

# n_step: number of steps in the cycle
# F_A, F_B: Free energy of A, B
# the amount of F particles consumed at each step (given by the CRN)

# Compare to steady state driving (flux of ATP)

using SymPy
using LinearAlgebra
using DifferentialEquations
using Catalyst
include("reactionNetworks.jl")
 
# Reaction: nA <-> B, where B is the higher free energy molecule
# Oscillate the volume between ϵ and 1/ϵ, where ϵ ≈ 1
# k_ratio = k2 / k1
# returns: N_B/N_A for the chosen parameters

### Equivalent ways of defining the ODE System # Going directly from the ReactionSystem: convert(ODESystem, rn)
# Using the S-matrix and j-space: f(du, u, p, t): du = S * j_space(N, rn, rates)  
# When using the latter method, need to divide certain columns to account for combinatoric rate law

# function optimalLandscape(nStates::Int64, α_max::Float64, A::Float64; inverting=true)
#     fwd = (inverting ? ((nStates-1)*α_max - A)/nStates : α_max)
#     rev = (inverting ? α_max : ((nStates-1)*α_max + A)/nStates)
#     landscape = repeat([fwd, rev], nStates)
#     inverting ? (landscape[end] = 0) : (landscape[end-1] = 0)
#     landscape
# end

function ratesFromDict(rn, rate_consts)
    pm = paramsmap(rn)
    rate_consts = sort(rate_consts, by=x->pm[x])
    rate_consts = collect(values(rate_consts))
    rate_consts
end

function mechanisticLandscape(additions, ΔG_AF, ΔG_FF, ΔG_A, intrinsic_barrier)
    steps = length(additions)
    landscape = zeros(Int64, steps * 2)
    landscape .= intrinsic_barrier

    for i in 1:steps-1
        landscape[2*i - 1] += additions[i]*ΔG_AF + (additions[i] - 1)*ΔG_FF
        i != 1 && (landscape[2*i-1] += ΔG_FF)
    end
    
    rem = ΔG_A - sum([landscape[2*i-1] - landscape[2*i] for i in 1:steps])
    (rem > 0) ? (landscape[end-1] += rem) : (landscape[end] -= rem)
    landscape
end


# Take in: reaction network, AF interaction, FF interaction, 4F->A free energy, inverse temperature
# Return: Dictionary of the reaction parameters
function setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, β, feed_rates; steps=4, intrinsic_barrier=3, r_0 = 100., comp = false, ΔG_BF = ΔG_AF, ΔG_B = ΔG_A)
    
    if (length(feed_rates) != length(reactionparams(rn)) - (comp+1)*2*steps)
        error("Length of feed rate vector is incorrect")
    end
    
    additions = [length(reactions(rn)[2*i-1].substoich) > 1 ? reactions(rn)[2*i-1].substoich[2] : 0 for i in 1:steps]
    comp && (additions_B = [length(reactions(rn)[2*i-1].substoich) > 1 ? reactions(rn)[2*i-1].substoich[2] : 0 for i in steps+1:2*steps])

    landscape = mechanisticLandscape(additions, ΔG_AF, ΔG_FF, ΔG_A, intrinsic_barrier)
    comp && (landscape = vcat(landscape, mechanisticLandscape(additions_B, ΔG_BF, ΔG_FF, ΔG_B, intrinsic_barrier)))

    k = r_0*exp.(-landscape*β)
    push!(k, feed_rates...)
    return landscape, Dict(zip(reactionparams(rn), k))
end

function singleReaction(k_ratio, stoich_a, A_tot, λ) 
    n = stoich_a
    S_h = [n*exp(-(n-1)*λ) n; exp(-(n-1)*λ) 1]
    S_l = [n*exp((n-1)*λ) n; exp((n-1)*λ) 1]
    S_avg = (S_h + S_l) / 2

    h_ss = nullspace(S_h)
    l_ss = nullspace(S_l)
    avg_ss = nullspace(S_avg)

    N_A = [solveset(Eq(n/k_ratio * exp(-(n-1)*α) * x^n + x, A_tot), x),
           solveset(Eq(n/k_ratio * exp((n-1)*α) * x^n + x, A_tot), x),
           solveset(Eq(n/k_ratio * cosh((n-1)*α) * x^n + x, A_tot), x)]
end

# Return: parameters of diff eqs, steady state solutions
function volumeOscillation(rn, λ0, Δλ, N0, rate_consts; t_end=10., additions=[])
    n = numreactions(rn); m = numspecies(rn)
    if length(rate_consts) != n || length(N0) != m
        error("Input parameters are of incorrect size")
    end
    
    rxs = reactions(rn)
    order = [sum(rxs[i].substoich) for i in 1:length(rxs)]
    order = [s <= 1 ? 0 : s-1 for s in order]
    if !isempty(additions)
        [order[i*2-1] = additions[i] for i in 1:length(additions)]
    end 

    if typeof(rate_consts) <: Dict 
        pm = paramsmap(rn)
        rate_consts = sort(rate_consts, by=x->pm[x])
        rate_consts = collect(values(rate_consts))
    end

    if typeof(N0) <: Dict 
        sm = speciesmap(rn)
        N0 = sort(N0, by=x->sm[x])
        N0 = collect(values(N0))
    end

    params_h = exp.(-order * (λ0+Δλ)) .* rate_consts
    params_l = exp.(-order * (λ0-Δλ)) .* rate_consts
    params_a = (params_l + params_h) / 2

    tspan = (0., t_end)
    rn_eq = convert(NonlinearSystem, rn, combinatoric_ratelaws = false)
    oprob_a = NonlinearProblem(rn_eq, N0, params_a); sol_a = solve(oprob_a, DynamicSS(Rodas5()))
    oprob_l = NonlinearProblem(rn_eq, N0, params_l); sol_l = solve(oprob_l, DynamicSS(Rodas5()))
    oprob_h = NonlinearProblem(rn_eq, N0, params_h); sol_h = solve(oprob_h, DynamicSS(Rodas5()))

    [params_l, params_h, params_a], [sol_l, sol_h, sol_a]
end

function j_space(N, rn, rates)
    n = numreactions(rn); J = zeros(n)
    rxs = reactions(rn); sm = speciesmap(rn)

    for i in 1:n
        reactants = rxs[i].substrates
        stoich = rxs[i].substoich
        idxs = [sm[r] for r in reactants]
        l = length(reactants)
        J[i] = rates[i] * prod([N[idxs[j]]^stoich[j] for j in 1:l])
    end

    J
end

function plotProd(sol_a, sol_l, sol_h, idx)
    a = [p[idx] for p in sol_a.u]
    l = [p[idx] for p in sol_l.u]
    h = [p[idx] for p in sol_h.u]
    plot([sol_a.t, sol_l.t, sol_h.t], [a, l, h], labels=["Average" "Low" "High"])
end

function sMatrix(rn, λ0, Δλ)
    n = numreactions(rn)
    m = numspecies(rn)
    S = netstoichmat(rn)
    S_l = float.(netstoichmat(rn))
    S_h = float.(netstoichmat(rn))

    for j in 1:n
        r_n = -sum(Int64[S[i,j] for i in 1:m if S[i,j] < 0])
        S_h[:,j] .*= exp(-(r_n - 1)*(λ0+Δλ))
        S_l[:,j] .*= exp(-(r_n - 1)*(λ0-Δλ))
    end

    S_avg = (S_l + S_h) / 2
    S_l, S_h, S_avg
end 

# By convention include the rate constants in this term. 
function A_V(rn, λ0, Δλ) 
    λ1 = λ0 + Δλ; λ2 = λ0 - Δλ
    rxs = reactions(rn)
    α = [sum(rx.substoich) for rx in rxs]

    A_eff = 0

    for i in 1:numreactions(rn)
        if i % 2 == 0 
            A_eff -= log(exp(-(α[i]-1)*λ1) + exp(-(α[i]-1)*λ2)/2)
        else
            A_eff += log(exp(-(α[i]-1)*λ1) + exp(-(α[i]-1)*λ2)/2)
        end
    end

    A_eff
end

function A_C(rn, ss, k_dict)
    A_c = 0
    rxs = reactions(rn)
    n = numreactions(rn)

    rates = zeros(n)
    pm = paramsmap(rn)
    for (k, r) in k_dict 
        rates[pm[k]] = r
    end
    sm = speciesmap(rn)

    for i in 1:numreactions(rn)
        rx = rxs[i]
        reactants = rx.substrates
        stoich = rx.substoich
        idxs = [sm[r] for r in reactants]
        if i % 2 == 0
            n_prod = rates[i]*prod([ss[idxs[j]]^stoich[j] for j in 1:length(idxs)])
            A_c -= log(n_prod) 
        else
            n_prod = rates[i]*prod([ss[idxs[j]]^stoich[j] for j in 1:length(idxs)])
            A_c += log(n_prod) 
        end 
    end

    A_c
end

# The affinity of the network, not the affinity of an individual autocatalyst undergoing a cycle. 
function A(rn, solutions, λ0, Δλ, k_dict)
    a_ss, l_ss, h_ss = [sol[end] for sol in solutions]

    A_va = A_V(rn, λ0, Δλ); A_ca = A_C(rn, a_ss, k_dict)
    A_vl = A_V(rn, λ0-Δλ, 0.); A_cl = A_C(rn, l_ss, k_dict)
    A_vh = A_V(rn, λ0+Δλ, 0.); A_ch = A_C(rn, h_ss, k_dict)
    return A_va+A_ca, A_vl+A_cl, A_vh+A_ch
end

# The affinity of an individual catalytic cycle, equal to A = -λΔσ - βΔF + \log (n_F^4)/n_A
function A(solutions, λ0, Δλ, ΔF; β = 1., Δσ = 3, m_idx = 1)
    l_ss, h_ss, a_ss = solutions
    A_a = Δσ*log((exp(-λ0-Δλ)+exp(-λ0+Δλ))/2) - β*ΔF + log(a_ss[2]^4 / a_ss[m_idx])
    A_l = -Δσ*(λ0-Δλ) - β*ΔF + log(l_ss[2]^4 / l_ss[m_idx])
    A_h = -Δσ*(λ0+Δλ) - β*ΔF + log(h_ss[2]^4 / h_ss[m_idx])
    A_l, A_h, A_a
end

function simulate(rn, ΔG_AF, ΔG_FF, ΔG_A, β, feed_rates, λ0, Δλ, N0; t_end = 10., comp = false)
    landscape, k = setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, 1., feed_rates, comp = comp)
    println(k)
    params, sols = volumeOscillation(rn, λ0, Δλ, N0, k);
    println(sols[1].u)
    aff = A(sols, λ0, Δλ, ΔG_A)
    sols, aff
end

function A_inst(rn, n_A, n_F, λ0, Δλ, ΔF; β=1.)
    S = sMatrix(rn, λ0, Δλ)
    affs = Float64[]
    for s in S
        V = log(s[1,1]*s[3,3]*s[4,5]*s[5,7] / (s[1,2]*s[3,4]*s[4,6]*s[5,8]))
        push!(affs, V-β*ΔF+log(n_F^4 / n_A))
    end
    affs
end

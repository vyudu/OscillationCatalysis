using LinearAlgebra
using Graphs, MetaGraphs
include("CycleOscillation.jl")

# Vertex properties
#   :energy
# Edge properties
#   :distance
#   :barrier
#   :ΔG (for electron reservoirs)
# barriers: Dict[Pair(Int64, Int64) => Float64]

function marcusRate(β::Float64, K0::Float64, ΔG::Float64, λ::Float64, d::Float64, decay=1.0)
    hbar = 6.5821 * 10^(-16)
    return 2*π/hbar * K0^2 * exp(-decay*d) * 1/(sqrt(4*π*λ/β))*exp(-β*(ΔG + λ)^2/4λ)
end

function singleMoleculeGraph(energies::Vector{<:Real}, barriers::Dict{Tuple{Int64, Int64},<:Real})
    nStates = length(energies)
    g = SimpleGraph(nStates)
    mg = MetaGraph(g)

    for edge in keys(barriers)
        add_edge!(mg, edge...)
        set_prop!(mg, edge..., :barrier, barriers[edge])
    end

    for i in 1:nStates
        set_prop!(mg, i, :energy, energies[i])
    end

    return mg
end

struct Cofactor
    name::String
    capacity::Int64
    potentials::Vector{<:Real}
end

# if the transfer goes to a resrvoir state, set the tuple as (idx, 0) where idx is the donor
function marcusGraph(cofactors::Vector{Cofactor}, distances::Dict{Tuple{Int64,Int64},<:Real}, resΔG::Dict{Int64, <:Real})
    nStates = foldr(*, [cf.capacity for cf in cofactors] .+ 1)
    capacities = [cf.capacity for cf in cofactors]
    potentials = [cf.potentials for cf in cofactors]
    nc = length(cofactors)
    g = SimpleGraph(nStates)
    mg = MetaGraph(g)

    # Add the transitions
    for i in 1:nStates
        st = state(i, nc)

        for (src, tgt) in keys(distances)
            dst = copy(st)
            
            if tgt == 0
                dst[src] -= capacities[src]
            else
                dst[src] -= 1; dst[tgt] += 1
            end

            if allowed(dst, capacities) 
                add_edge!(mg, i, index(dst, nc))
                set_prop!(mg, i, index(dst, nc), :distance, distances[(src, tgt)])
                tgt == 0 && set_prop!(mg, i, index(dst, nc), :ΔG, resΔG[src])
            end
        end
    end
    
    # Add the energies and the states
    for i in 1:nStates
        st = state(i, nc); energy = 0
        energy = sum([sum(potentials[i][1:st[i]]) for i in 1:nc])
        set_prop!(mg, i, :energy, energy)
        set_prop!(mg, i, :state, state(i, nc))
    end

    return mg
end

# Helpers to convert from the vertex index to the state
function state(idx::Int64, nc::Int64)
    state = zeros(Int, nc) 
    idx -= 1

    for i in 1:nc
        state[i] = floor(Int, idx / (2^(nc-i)))
        idx -= 2^(nc-i)*state[i]
    end
    return state
end

function index(state::Vector{Int64}, nc::Int64)
    return sum([2^(nc-i) * state[i] for i in 1:nc]) + 1
end

function allowed(state::Vector{Int64}, capacities::Vector{Int64})
    return all(0 .<= state .<= capacities)
end

function rateMatrix(G::MetaGraph, β::Float64; marcus=false, λ=1.0, K0=1.0) 
    R = zeros(nv(G), nv(G))

    for v in vertices(G)
        for n in neighbors(G, v)
            if marcus
                ΔG = 0 
                if has_prop(G, v, n, :ΔG)
                    v < n ? ΔG = get_prop(G, v, n, :ΔG) : ΔG = -get_prop(G, v, n, :ΔG) 
                else
                    ΔG = get_prop(G, n, :energy) - get_prop(G, v, :energy)
                end
                d = get_prop(G, v, n, :distance)
                R[n, v] = marcusRate(β, K0, ΔG, λ, d)
            else
                B = get_prop(G, v, n, :barrier)
                E = get_prop(G, v, :energy)
                R[n, v] = rate(β, E, B) 
            end
        end
        totalRate = sum(R[:, v])
        R[v, v] = -totalRate
    end

    return R
end

function pop(P::Vector{Float64}, idx::Int64, nc::Int64, electrons::Int64) 
    tot = 0
    for i in 1:length(P) 
        st = state(i, nc)
        (st[idx] == electrons) && (tot += P[i]) 
    end
    return tot 
end

function flux(R::Matrix{Float64}, P::Vector{Float64}, cofactors::Vector{Cofactor}, idx::Int64)
    nc = length(cofactors)
    ox_i = cofactors[idx].capacity
    ox_f = 0
    st_i = zeros(Int, nc); st_i[idx] = ox_i
    st_f = zeros(Int, nc); st_f[idx] = ox_f
    i = index(st_i, nc); f = index(st_f, nc) 

    return (pop(P, idx, nc, ox_i)*R[f,i] - pop(P, idx, nc, ox_f)*R[i,f]) * ox_i
end

function fluxPlot(cofactors::Vector{Cofactor}, distances::Dict{Tuple{Int64,Int64},<:Real}, resΔG::Dict{Int64,<:Real}, ΔG_range::AbstractRange; β=38.91, K0=0.1, T=10.) 
    fluxD = Float64[]; fluxH = Float64[]; fluxL = Float64[];
    for ΔG in ΔG_range
        resΔG[1] = ΔG
        mg = marcusGraph(cofactors, distances, resΔG)
        R = rateMatrix(mg, β, marcus=true, K0=K0)
        P_i = zeros(size(R)[1]); P_i[1] = 1.

        (d,h,l) = sort(collect(keys(resΔG)))
        push!(fluxD, flux(R, evolve(P_i, T, R), cofactors, d))
        push!(fluxH, flux(R, evolve(P_i, T, R), cofactors, h))
        push!(fluxL, flux(R, evolve(P_i, T, R), cofactors, l))
    end
    plt = plot(ΔG_range, [fluxD, fluxH, fluxL])
end

function cofactorArray(branchLength::Int64, slope::Float64, E_B::Float64)
    B = Cofactor("B", 2, [-E_B, E_B])
    H = Cofactor[]; L = Cofactor[];
    for i in 1:branchLength
        push!(H, Cofactor("H$i", 1, [E_B - slope*i]))
        push!(L, Cofactor("L$i", 1, [-E_B + slope*i]))
    end
    return [B, H..., L...]
end

function distanceOscillation(cofactors::Vector{Float64}, distances::Dict{Tuple{Int64,Int64},<:Real}, resΔG::Dict{Int64,<:Real}, dilation::Float64; β=38.91, K0=0.1, T=10.)
    d_i = Dict(edge => dist/dilation for (edge,dist) in distances)
    d_f = Dict(edge => dist*dilation for (edge, dist) in distances) 

    mg_i = marcusGraph(cofactors, d_i, resΔG)
    mg_f = marcusGraph(cofactors, d_f, resΔG)
    mg_ss = marcusGraph(cofactors, distances, resΔG)

    R_i = rateMatrix(mg_i, β, marcus=true, K0=K0)
    R_f = rateMatrix(mg_f, β, marcus=true, K0=K0)
    R_avg = (R_i + R_f) / 2
    R_ss = rateMatrix(mg_ss, β, marcus=true, K0=K0)
    P_i = zeros(size(R_i)[1]); P_i[1] = 1.

    flux_i = flux(R_i, evolve(P_i, T, R_i), cofactors, 1)
    flux_f = flux(R_f, evolve(P_i, T, R_f), cofactors, 1)
    flux_avg = flux(R_avg, evolve(P_i, T, R_avg), cofactors, 1)
    flux_ss = flux(R_ss, evolve(P_i, T, R_ss), cofactors, 1)
    Dict("Stationary"=>flux_ss, "Oscillating"=>flux_avg, "Low"=>flux_f, "High"=>flux_i)
end

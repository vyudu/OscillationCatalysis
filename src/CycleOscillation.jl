# We look at cycles, where rate matrices can be expressed as vectors
# Landscapes: [E1, E2, ..., En, B12, B23, ..., Bn1]
# Rate vectors: [r21, r32, ..., r1n, r12, r23, ..., rn1]

# ===== DYNAMICS ===== #

function rate(β, E, B)
    return exp(-β*(B - E))
end

function optimalLandscape(nStates, α_max, A; inverting=true)
    fwd = (inverting ? ((nStates-1)*α_max - A)/nStates : α_max)
    rev = (inverting ? α_max : ((nStates-1)*α_max + A)/nStates)
    Δ = rev - fwd

    energies = [n*Δ for n in (nStates-1):-1:0]
    barriers = [E + fwd for E in energies]
    inverting ? (barriers[end] = energies[1]) : (barriers[end] = energies[end])
    return [energies..., barriers...]
end

function buildRateMatrix(β, landscape::Array{T}, A) where T<: Real
    nStates = Int(length(landscape) / 2)
    rateMatrix = zeros(nStates, nStates)
    energies = landscape[1 : nStates]
    barriers = landscape[nStates+1 : 2*nStates]
        
    for i in 1:nStates-1
        rateMatrix[i, i+1] = exp(-β*(barriers[i] - energies[i+1]))
        rateMatrix[i+1, i] = exp(-β*(barriers[i] - energies[i]))
    end
    
    rateMatrix[1, nStates] = exp(-β*(barriers[nStates] - energies[nStates] - A))
    rateMatrix[nStates, 1] = exp(-β*(barriers[nStates] - energies[1]))

    for i in 1:nStates
        totalRate = sum(rateMatrix[:, i])
        rateMatrix[i, i] = -totalRate
    end
    
    return rateMatrix 
end

function evolve(previous, timestep, R)
    return exp(R*timestep) * previous
end

# ===== CYCLE PROPERTIES ===== #

function affinity(R::Matrix{T}) where T <: Real
    nStates = size(R)[1]
    return sum( [ log(R[i+1, i] / R[i, i+1]) for i in 1:nStates-1]) + log(R[1, nStates] / R[nStates, 1])
end

function affinity(r::Array{T}) where T <: Real
    nStates = Int(size(r)[1] / 2)
    return sum( [log(r[i]/r[i + nStates]) for i in 1:nStates])
end

function matrixToVector(R, nStates)
    p = [R[i+1, i] for i in 1:nStates-1]
    push!(p, R[1, nStates])
    m = [R[i, i+1] for i in 1:nStates-1]
    push!(m, R[nStates, 1])
    return vcat(p, m)
end

function current(r::Array{T}, nStates::Int64) where T<:Real
    num = foldr(*, r[1:nStates]) - foldr(*, r[nStates+1:2*nStates])
    denom = sum([matrixTree(r, nStates, n) for n in 1:nStates])
    return num / denom
end

function matrixTree(r::Array{T}, nStates::Int64, n) where T<: Real
    denom = 0
    for i in 1:nStates
        if i == n
            term = foldr(*, r[nStates+1:2*nStates]) / r[i == 1 ? nStates + nStates : nStates + i - 1]
            denom += term
        elseif i > n
            term = foldr(*, r[i:nStates]) * foldr(*, r[1:n - 1]) * foldr(*, r[nStates+n : nStates+i-2])
            denom += term
        elseif i < n
            term = foldr(*, r[nStates+1 : nStates+i-2]) * foldr(*, r[i:n-1]) * foldr(*, r[nStates+n:2*nStates])
            denom += term
        end
    end
    return denom
end

function curvature(α)
    nStates = len(α) / 2
    curvature = foldr(*, α[1:nStates].^2) - foldr(*, α[nStates+1:2*nStates].^2)
    return curvature
end


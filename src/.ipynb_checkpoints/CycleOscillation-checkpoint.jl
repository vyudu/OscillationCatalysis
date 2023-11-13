using LinearAlgebra
using QuadGK

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


# ===== PLOTS ===== #


function tempOscillation(landscape, signal, β1, β2, A, nStates)
    R1 = buildRateMatrix(β1, landscape, A, nStates)
    R2 = buildRateMatrix(β2, landscape, A, nStates)
    
    # p = eigvecs(R1)[:,nStates] 
    p = rand(nStates)
    p = p ./ sum(p) 

    currAvgs = Float64[]

    c1(p, t) = R1[1,nStates]*(exp(R1*t)*p)[nStates] - R1[nStates,1]*(exp(R1*t)*p)[1]
    c2(p, t) = R2[1,nStates]*(exp(R2*t)*p)[nStates] - R2[nStates,1]*(exp(R2*t)*p)[1]

    for int in signal
        curr1, err = quadgk(t->c1(p,t), 0, int/2, rtol=1e-8)
        p = evolve(p, int/2, R1)
        curr2, err = quadgk(t->c2(p,t), 0, int/2, rtol=1e-8)
        p = evolve(p, int/2, R2)

        push!(currAvgs, ((curr1+curr2)/int))
    end

    transitionTimes = accumulate(+, signal)
    pushfirst!(transitionTimes, 0)
    β(t) = any([0 <= (t-transitionTimes[i]) < (signal[i]/2) for i in 1:length(signal)]) ? β1 : β2
    
    p1 = plot(β, xlims=(0,transitionTimes[end]), title="Temperature Signal", xlabel="time", ylabel="β", fontsize=10)
    p2 = scatter(transitionTimes[2:end], currAvgs, title="Period-averaged Current", xlabel="time", ylabel="flux", fontsize=10)
    return p1, p2
end


function steadyStateCurrent(period, R1, R2; tol=1e-8)
    nStates = size(R1)[1]
    W = exp(R2*period/2)*exp(R1*period/2)
    p = eigvecs(W-I)[:,3] ./ sum(eigvecs(W-I)[:,3])

    # Integrate the current over a full period
    c1(p, t) = R1[1,nStates]*(exp(R1*t)*p)[nStates] - R1[nStates,1]*(exp(R1*t)*p)[1]
    c2(p, t) = R2[1,nStates]*(exp(R2*t)*p)[nStates] - R2[nStates,1]*(exp(R2*t)*p)[1]

    int1, err = quadgk(t->c1(p,t), 0, period/2, rtol=1e-8)
    p = exp(R1*period/2)*p
    int2, err = quadgk(t->c2(p,t), 0, period/2, rtol=1e-8)
    
    return (int1 + int2) / period
end


function currentVsFrequency(landscape, β1, β2, F_B, range::AbstractRange, nStates)
    R1 = buildRateMatrix(β1, landscape, F_B, nStates)
    R2 = buildRateMatrix(β2, landscape, F_B, nStates)
    freqs = Float64[]
    currents = Float64[]

    for period in exp.(range)
        curr = steadyStateCurrent(period, R1, R2)
        push!(currents, curr)
        push!(freqs, 1/period)
    end

    return freqs, currents
end

function ΔAvsCurvature()

end

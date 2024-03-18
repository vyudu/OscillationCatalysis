using OscillationCatalysis, QuadGK, LinearAlgebra, Plots

β0 = 1.0; Δβ = 0.5; ΔF = 1.0
α_max = 8.; nStates = 3 
inverting = optimalLandscape(nStates, α_max, 1.)
assisting = optimalLandscape(nStates, α_max, 1., inverting=false)

### AFFINITY VS CURVATURE PLOTS ###

function generateAlpha(scale, aff=1)
    alphas = rand(5) * scale

    while !(0 < sum(alphas[1:3]) - alphas[4] - alphas[5] + aff < scale )
        alphas = rand(5) * scale
    end
    return [alphas..., sum(alphas[1:3]) - alphas[4] - alphas[5] + aff]
end

# a21, a32, a13, a12, a23, a31
function generateLandscapeFromAlpha(alphas)
    E1 = 0
    B12 = alphas[1]
    E2 = B12 - alphas[4]
    B23 = E2 + alphas[2]
    E3 = B23 - alphas[5]
    B31 = alphas[6]

    return [E1, E2, E3, B12, B23, B31]
end

function higherOrderTerm(δβ, alphas, order)
    return δβ^2 * (sum(alphas[1:3] .^ order) - sum(alphas[4:6] .^ order)) / 12
end

# Given a pair of temperatures, free energy difference, and number of points, randomly generate landscapes and compute their affinity and curvatures
function linearApproximation(β1, β2, ΔF, scale, points)
    Δaff = []
    curvatures = []
    higherOrderTerms = []
    β0 = (β1 + β2) / 2

    for i in 1:points
        α = generateAlpha(scale)
        landscape = generateLandscapeFromAlpha(α)

        R1 = buildRateMatrix(β1, landscape, ΔF, 0)
        R2 = buildRateMatrix(β2, landscape, ΔF, 0)
        R_avg = (R1 + R2) / 2

        aff = affinity(R_avg)

        push!(Δaff, aff - ΔF * β0)
        push!(curvatures, curvature(α))
        push!(higherOrderTerms, higherOrderTerm((β1 - β2) / 2, α, 4))

    end

    return Δaff, curvatures, higherOrderTerms
end

function ΔAvsCurvature(β0, Δβ, ΔF, scale, points) 
    aff, curv, fot = linearApproximation(β0-Δβ, β0+Δβ, ΔF, scale, points)
    s1 = scatter(curv, aff/(Δβ^2)*2, xlabel="Curvature", ylabel="2ΔA / Δβ^2", xtickfontsize=10, ytickfontsize=10, legend=false, aspect_ratio=1)
end



### CURRENT VS TIME PLOT ###

function tempOscillation(signal; nStates=3)
    R1 = buildRateMatrix(β0-Δβ, inverting, ΔF) 
    R2 = buildRateMatrix(β0+Δβ, inverting, ΔF)
    
    p = eigvecs(R2)[:,nStates] 
    p = p ./ sum(p) 

    currAvgs = Float64[]

    for period in signal
        curr, p = cycleCurrent(period, R1, R2; p = p)
        push!(currAvgs, curr)
    end

    transitionTimes = accumulate(+, signal)
    pushfirst!(transitionTimes, 0)
    β(t) = any([0 <= (t-transitionTimes[i]) < (signal[i]/2) for i in 1:length(signal)]) ? β0-Δβ : β0+Δβ
    samplePts = 0:(signal[end]/100):transitionTimes[end]
    
    p1 = plot(samplePts, β, xlims=(0,transitionTimes[end]), title="Temperature Signal", xlabel="time", ylabel="β", fontsize=10)
    p2 = scatter(transitionTimes[2:end], currAvgs, xlims=(0, transitionTimes[end]), title="Period-averaged Current", xlabel="time", ylabel="flux", fontsize=10)
    return p1, p2
end


### CURRENT VS FREQUENCY PLOTS ###

function cycleCurrent(period, R1, R2; tol=1e-8, p=zeros(nStates))
    nStates = size(R1)[1]

    # Assume steady state if initial p is not provided
    if p == zeros(nStates)
        W = exp(R2*period/2)*exp(R1*period/2)
        p = eigvecs(W-I)[:,nStates] ./ sum(eigvecs(W-I)[:,nStates])
    end

    # Integrate the current over a full period
    c1(p, t) = R1[1,nStates]*(exp(R1*t)*p)[nStates] - R1[nStates,1]*(exp(R1*t)*p)[1]
    c2(p, t) = R2[1,nStates]*(exp(R2*t)*p)[nStates] - R2[nStates,1]*(exp(R2*t)*p)[1]

    curr1, err = quadgk(t->c1(p,t), 0, period/2, rtol=1e-8)
    p = exp(R1*period/2)*p
    curr2, err = quadgk(t->c2(p,t), 0, period/2, rtol=1e-8)
    p = exp(R2*period/2)*p
    
    return (curr1+curr2) / period, p
end

function currentVsFrequency(range::AbstractRange; nStates=3)
    R1i = buildRateMatrix(β0-Δβ, inverting, ΔF)
    R2i = buildRateMatrix(β0+Δβ, inverting, ΔF)
    R1a = buildRateMatrix(β0-Δβ, assisting, ΔF)
    R2a = buildRateMatrix(β0+Δβ, assisting, ΔF)

    freqs = Float64[]
    currents_a = Float64[]
    currents_i = Float64[]

    for freq in exp.(range)
        curr_a, p = cycleCurrent(1/freq, R1a, R2a)
        curr_i, p = cycleCurrent(1/freq, R1i, R2i)
        push!(currents_a, curr_a)
        push!(currents_i, curr_i)
        push!(freqs, freq)
    end

    plot(log.(freqs), [currents_a, currents_i], xlabel="log(f)", ylabel="Current", labels=["Assisting" "Inverting"])
    title!("Reaction Flux vs. Oscillation Frequency for optimal landscapes")
end

include("CycleOscillation.jl")

α_max = 10.; A = 1., nStates = 3
inverting = optimalLandscape(α_max, nStates, A)
assisting = optimalLandscape(α_max, nStates, A)

freqs_i, currents_i = currentVsFrequency(inverting, β1, β2, A, -4:0.01:2, nStates)
freqs_a, currents_a = currentVsFrequency(assisting, β1, β2, A, -4:0.01:2, nStates)

fig1 = plot(freqs_i, [currents_i, currents_a])
savefig(fig1, "PRE22fig1.jpg")

signal = [1. * 0.99^i for i in 1:800]
fig2, fig3 = tempOscillation(inverting, signal, β1, β2, A, nStates)

#TODO: fig 4 - curvature vs affinity diff


include("CRNOscillation.jl")
include("reactionNetworks.jl")

# TODO: 
### Affinity vs. lambda for a single network
### Affinity vs. Curly C for a single network (smaller oscillations)

# Define Parameters
# ΔG_A is the energy associated with reaction 4F → A
# ΔG_AF is the interaction energy between an A and an F
ΔG_AF = 1; ΔG_BF = 1; ΔG_A = 6; ΔG_B = 6; ΔG_FF = 1
β = 1.; λ0 = 0; Δλ = 1.5 

# Species: A, F, AF, AF2, AF3, B, (BF2/BB), (BF3/BB/BBa), (BF4/BB'/BBb)
N0 = [.1, 1., 0., 0., 0., 0., .1, 0.]

# Energy landscape: k_1, k_2, k_3, k_4
landscape_4 = mechanisticLandscape([1,1,1,1], ΔG_AF, ΔG_FF, ΔG_A) 
landscape_3 = mechanisticLandscape([2,1,1,0], ΔG_AF, ΔG_FF, ΔG_A)
landscape_2 = mechanisticLandscape([2,2,0,0], ΔG_AF, ΔG_FF, ΔG_A)
landscape_1 = mechanisticLandscape([4,0,0,0], ΔG_AF, ΔG_FF, ΔG_A)
# Rates: feed reactions
r_f = [10., 5., 10., 5., 10., 5.]
r_0 = 100.
rates_c = vcat(r_0*exp.(-landscape_4), r_0*exp.(-landscape_1))
rates_o = vcat(rates_c, r_f) 

# Plot the change in the exit rate of A, B as a function of time, turning on oscillation halfway through
function competitionPlot() 
    params_l, params_h, params_a = volumeOscillation(competitive, λ0, Δλ, N0, rates_o)
    ode = convert(ODESystem, competitive)
    tspan = (0., 500.)
    prob = ODEProblem(ode, N0, tspan, params_h)
    prob2 = ODEProblem(ode, N0, tspan, params_l)

    # begin oscillation at time t = 50.
    affect!(integrator) = integrator.p = params_a
    cb = PresetTimeCallback([tspan[2] / 2], affect!)
    sol = solve(prob, Tsit5(), callback = cb)
    sol2 = solve(prob2, Tsit5(), callback = cb)

    lr_A = [p[1] for p in sol.u] 
    lr_B = [p[6] for p in sol.u]
    lr_A_l = [p[1] for p in sol2.u] 
    lr_B_l = [p[6] for p in sol2.u]

    p = plot(sol.t, [lr_A,lr_B], labels = ["Leaving rate of A" "Leaving rate of B"])
    plot!(p, sol2.t, [lr_A_l, lr_B_l], labels = ["Leaving rate of A (high)" "Leaving rate of B (high"])
    sol, sol2, p
end

function competition_closed()
    params_l, params_h, params_a = volumeOscillation(competitive_closed, λ0, Δλ, N0, rates_c)
    ns = convert(NonlinearSystem, competitive_closed)
    ss_prob = NonlinearProblem(ns, N0, params_l); ss = solve(ss_prob, DynamicSS(Rodas5()))
    ss_prob_h = NonlinearProblem(ns, N0, params_h); ss_h = solve(ss_prob_h, DynamicSS(Rodas5()))

    ode = convert(ODESystem, competitive_closed)
    tspan = (0., 50.)
    affect!(integrator) = integrator.p = params_a
    cb = PresetTimeCallback([tspan[2] / 2], affect!)
    prob = ODEProblem(ode, ss.u, tspan, params_l); sol = solve(prob, Tsit5(), callback=cb)
    prob_h = ODEProblem(ode, ss_h.u, tspan, params_h); sol_h = solve(prob_h, Tsit5(), callback=cb)

    A = [p[1] for p in sol.u]
    B = [p[6] for p in sol.u]
    A_h = [p[1] for p in sol_h.u]
    B_h = [p[6] for p in sol_h.u]
    p = plot([sol.t], [A,B], labels = ["A concentration" "B concentration"]) 
    plot!(p, [sol_h.t], [A_h,B_h], labels = ["A concentration (high)" "B concentration (high)"]) 
    sol, sol_h, p
end

function affinityPlot(open=true)
    landscapes = [landscape_1, landscape_2, landscape_3, landscape_4]; rns = [oneStep_closed, twoStep_closed, threeStep_closed, autocatalysis_closed];
    rates = [setRateConstants(landscape, rn, 1.) for (landscape, rn) in zip(landscapes, rns)]
    params4, sols4 = volumeOscillation(autocatalysis_closed, λ0, Δλ, N0[1:5], rates[4])
    params3, sols3 = volumeOscillation(threeStep_closed, λ0, Δλ, N0[1:5], rates[3]) 
    params2, sols2 = volumeOscillation(twoStep_closed, λ0, Δλ, N0[1:5], rates[2])
    params1, sols1 = volumeOscillation(oneStep_closed, λ0, Δλ, N0[1:5], rates[1])

    ΔA = zeros(Float64, 0)
    for sol in [sols1, sols2, sols3, sols4]
        A_l, A_h, A_a = A(sol, λ0, Δλ, ΔG_A)
        push!(ΔA, A_a - A_l)
    end

    p = bar(collect(1:4), ΔA, title="Affinity change upon oscillation", xlabel="Steps in CRN", ylabel="ΔA")  
    savefig(p, "~/workspace/OscillationCatalysis/plots/affinity.png")
end

function pumpedPlot(num="")
    params_l, params_h, params_a = volumeOscillation(autocatalysis_open, λ0, Δλ, N0[1:5], [rates_o[1:8]..., 1., 1.])
    ode = convert(ODESystem, competitive)
    tspan = (0., 50.)
    prob = ODEProblem(ode, N0, tspan, params_l)

    # begin oscillation at time t = 50.
    affect!(integrator) = integrator.p = params_a
    cb = PresetTimeCallback([tspan[2] / 2], affect!)
    sol = solve(prob, Tsit5(), callback = cb)

    lr_A = [p[1] for p in sol.u] 

    p = plot(sol.t, [lr_A], labels = ["Leaving rate of A"])
    savefig(p, "BioCosmosPlots/pumpedPlot"*num*".png")
    sol, p
end

function frequencyPlots()
    affect!(integrator) = (integrator.p == params_l) ? params_h : params_l
    cb = PresetTimeCallback([tspan[2]])
end

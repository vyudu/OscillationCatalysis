using OscillationCatalysis, Catalyst, DifferentialEquations, Sundials
const o = OscillationCatalysis

# Define Parameters
# ΔG_A is the energy associated with reaction 4F → A
# ΔG_AF is the interaction energy between an A and an F
ΔG_AF = 1; ΔG_BF = 1; ΔG_A = 6; ΔG_B = 6; ΔG_FF = 1
β = 1.; λ0 = 0; Δλ = 0.5 

# Species: A1, F, A2, A3, A4, B1, B2, B3, B4
N0 = [.1, 1., 0., 0., 0., .1, 0., 0., 0.]
init_c = Dict(zip(species(o.competitive_A4B1_open), N0))
init_nc = Dict(zip(species(o.autocatalysis_open), N0[1:5]))

# Rates: feed reactions
feed_rates_c = [8., 1., 1., 1., 1., 1.]
feed_rates_nc = [4., 1., 1., 1.]

# Networks
comp_rns = [o.competitive_A4B1_open, o.competitive_A4B2_open, o.competitive_A4B3_open]
single_rns = [o.oneStep_open, o.twoStep_open, o.threeStep_open, o.autocatalysis_open];

# Case 1: fix N_A, N_F, compute instantaneous affinity
function A_inst_bargraph(n_a, n_f) 
   l = A_inst(oneStep_open, n_a, n_f, λ0, Δλ, ΔG_A)[1]
   h = A_inst(oneStep_open, n_a, n_f, λ0, Δλ, ΔG_A)[2]

   affs = [A_inst(rn, n_a, n_f, λ0, Δλ, ΔG_A)[3] for rn in single_rns]  
   p = bar(collect(1:4), affs, legend=false)
   xlabel!("Steps in Reaction Network")
   ylabel!("Instantaneous Affinity")

   Plots.abline!(0, h, color=:blue); Plots.abline!(0, l, color=:red); Plots.abline!(0,0, color=:black)
   ylims!((floor(h), ceil(l)))
   title!("Instantaneous Affinity vs. Number of Steps")
   annotate!([(4, h-0.2, "High Volume", 8), (4, l-0.2, "Low Volume", 8)])
end

# Case 2: fix N_F, compute the "non-competitive" fitness (proportion that goes into A, B at equilibrium)
function nc_fitness(rn, n_l, n_h) 
    ratios_l = Float64[]
    ratios_a = Float64[]
    ratios_h = Float64[]

    ls, rate_consts = o.setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, β, feed_rates_c, comp=true)
    k = ratesFromDict(rn, rate_consts)
    additions = [length(reactions(rn)[2*i-1].substoich) > 1 ? reactions(rn)[2*i-1].substoich[2] : 0 for i in 1:8]

    ratio_A(sol) = sum(sol[1:4]) / sum(sol)
    ratio_B(sol) = sum(sol[5:8]) / sum(sol)

    for n_F in range(n_l, n_h, 100)
        k_primes = k[1:16]
        [k_primes[2*i-1] *= n_F^(additions[i]) for i in 1:8]

        N0 = [n_F / 4, 0., 0., 0., n_F /4, 0., 0., 0.]
        params, sols = o.volumeOscillation(rn, λ0, Δλ, N0, k_primes, additions=additions)
        push!(ratios_l, ratio_B(sols[1].u))
        push!(ratios_h, ratio_B(sols[2].u))
        push!(ratios_a, ratio_B(sols[3].u))
    end

    plot(range(n_l, n_h, 100), [ratios_l, ratios_h, ratios_a], labels=["Low Volume" "High Volume" "Average Volume"], title="Non-competitive Fitness", legendfontsize=16, titlefontsize=18, tickfontsize=14)
end

# Plot the concentrations of each species in the reactor as a function of time
function competition(rn; tspan = (0., 100.), switch = tspan[2] / 2) 
    landscape, k = setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, 1., feed_rates_c, comp = true)
    params, sols = volumeOscillation(rn, λ0, Δλ, init_c, k)
    
    ode = convert(ODESystem, rn, combinatoric_ratelaws = false)
    prob_l = ODEProblem(ode, N0, tspan, params[1]) # Start out at the small volume
    prob_h = ODEProblem(ode, N0, tspan, params[2])

    # begin oscillation at time t = 50.
    affect!(integrator) = integrator.p = params[3]
    cb = PresetTimeCallback(switch, affect!)
    sol_l = solve(prob_l, CVODE_BDF(), callback = cb)
    sol_h = solve(prob_h, CVODE_BDF(), callback = cb)

    lr_A = [p[1] for p in sol_l.u] 
    lr_B = [p[6] for p in sol_l.u]
    lr_A_h = [p[1] for p in sol_h.u] 
    lr_B_h = [p[6] for p in sol_h.u]

    p = plot(sol_l.t, [lr_A,lr_B], labels = ["Leaving rate of A" "Leaving rate of B"])
    plot!(sol_h.t, [lr_A_h, lr_B_h], labels = ["Leaving rate of A (high)" "Leaving rate of B (high)"])
    sol_l, sol_h, p
end

# Plot the steady state population of A, B as a function of time in a closed network 
function competition_closed(rn, file_name; tspan = (0., 100.))
    landscape, k = setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, 1., feed_rates_c, comp = true)
    params, sols = volumeOscillation(rn, λ0, Δλ, N0, k)

    ode = convert(ODESystem, competitive_closed)
    affect!(integrator) = integrator.p = params_a
    cb = PresetTimeCallback([tspan[2] / 2], affect!)
    
    prob = ODEProblem(ode, ss.u, tspan, params[1]); sol = solve(prob, Tsit5(), callback=cb)
    prob_h = ODEProblem(ode, ss_h.u, tspan, params[2]); sol_h = solve(prob_h, Tsit5(), callback=cb)

    A = [p[1] for p in sol.u]
    B = [p[6] for p in sol.u]
    A_h = [p[1] for p in sol_h.u]
    B_h = [p[6] for p in sol_h.u]
    p = plot([sol.t], [A,B], labels = ["A concentration" "B concentration"]) 
    plot!(p, [sol_h.t], [A_h,B_h], labels = ["A concentration (high)" "B concentration (high)"]) 
    sol, sol_h, p
end

# Plot the affinity of a given non-competitive cycle as the number of steps increases from one to four. 
function A_numSteps()
    rates = [setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, β, feed_rates_nc) for rn in single_rns]
    solutions = [volumeOscillation(rn, λ0, Δλ, init_nc, k) for (rn, (ls, k)) in zip(single_rns, rates)]
    avg_solutions = [volumeOscillation(rn, λ0, 0., init_nc, k) for (rn, (ls, k)) in zip(single_rns, rates)]
    
    ΔA_l = zeros(Float64, 0)
    ΔA_h = zeros(Float64, 0)
    ΔA_a = zeros(Float64,0)
    
    for ((params, sols), (avgparams, avgsols)) in zip(solutions, avg_solutions)
        A_l, A_h, A_a = A(sols, λ0, Δλ, ΔG_A)
        A_avg = A(avgsols, λ0, 0., ΔG_A)[1]
        println("low: $A_l high: $A_h avg: $(A_a-A_avg)")
        push!(ΔA_a, A_a - A_avg)
        # push!(ΔA_h, A_h - A_a) 
    end

    p = groupedbar(ΔA_a, title="Affinity change upon oscillation", xlabel="Steps in CRN", ylabel="ΔA")  
    p
end

function A_numSteps_comp()
    rates = [setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, β, feed_rates_c, comp=true) for rn in comp_rns]
    solutions = [volumeOscillation(rn, λ0, Δλ, init_c, k) for (rn, (ls, k)) in zip(comp_rns, rates)]
    
    ΔA_l = zeros(Float64, 0)
    ΔA_h = zeros(Float64, 0)
    ΔA_a = zeros(Float64, 0)
    
    for (params, sols) in solutions
        A_l, A_h, A_a = A(sols, λ0, Δλ, ΔG_A)
        A_lb, A_hb, A_ab = A(sols, λ0, Δλ, ΔG_A, m_idx=6)
        println("low: $A_l high: $A_h avg: $A_a")
        println("low: $A_l high: $A_h avg: $A_a")
        push!(ΔA_l, A_l - A_lb)
        push!(ΔA_h, A_h - A_hb) 
        push!(ΔA_a, A_a - A_ab) 
    end

    p = groupedbar([ΔA_l ΔA_h ΔA_a], label=["low" "high"], title="Affinity change upon oscillation", xlabel="Steps in CRN", ylabel="ΔA-ΔB")  
    p
end

function lr_numSteps()
    rns = [oneStep_open, twoStep_open, threeStep_open, autocatalysis_open];
    rates = [setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, β, feed_rates_nc) for rn in rns]
    solutions = [volumeOscillation(rn, λ0, Δλ, init_nc, k) for (rn, (ls, k)) in zip(rns, rates)]
    
    lrs_l = zeros(Float64, 0)
    lrs_h = zeros(Float64, 0)
    lrs_a = zeros(Float64, 0)
    
    for (params, sols) in solutions
        lrs = [sols[i].u[1] for i in 1:3]
        push!(lrs_l, lrs[1]); push!(lrs_h, lrs[2]); push!(lrs_a, lrs[3])
    end

    p = groupedbar([lrs_l lrs_h lrs_a], label=["low" "high" "avg"], title="Leaving rate change upon oscillation", xlabel="Steps in CRN", ylabel="Leaving Rate")  
    p
end

function lr_numSteps_comp()
    rns = [competitive_A4B3_open, competitive_A4B2_open, competitive_A4B1_open];
    rates = [setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, β, feed_rates_c, comp=true) for rn in rns]
    solutions = [volumeOscillation(rn, λ0, Δλ, init_c, k) for (rn, (ls, k)) in zip(rns, rates)]
    
    lrs_l = zeros(Float64, 0)
    lrs_h = zeros(Float64, 0)
    lrs_a = zeros(Float64, 0)
    
    for (params, sols) in solutions
        lrs_A = [sols[i].u[1] for i in 1:3]
        lrs_B = [sols[i].u[6] for i in 1:3]
        lrs = lrs_A .- lrs_B
        push!(lrs_l, lrs[1]); push!(lrs_h, lrs[2]); push!(lrs_a, lrs[3])
    end

    p = groupedbar([lrs_l lrs_h lrs_a], label=["low" "high" "avg"], title="Leaving rate change upon oscillation", xlabel="Steps in B's cycle", ylabel="LR_A - LR_B")  
    p
end

# Plot the affinity of a single non-competitive network as a function of the log-volume. 
function A_λ(rn, λ_l, λ_h)
    landscape, k = setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, 1., feed_rates_nc, comp = false)
    affs = Float64[]
    lrs = Float64[]
    
    for λ in range(λ_l, λ_h, 100)
        params, sols = volumeOscillation(rn, λ, 0., init_nc, k)
        lr = sols[1].u[1]
        aff = A(sols, λ, 0., ΔG_A)
        push!(lrs, lr)
        push!(affs, aff[1])
    end
    p = plot(collect(range(λ_l, λ_h, 100)), [affs], title = "Affinity vs. λ", labels = "Affinity")
    p 
end

function lr_λ(rn, λ_l, λ_h)
    landscape, k = setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, 1., feed_rates_nc, comp = false)
    affs = Float64[]
    lrs = Float64[]
    
    for λ in range(λ_l, λ_h, 100)
        params, sols = volumeOscillation(rn, λ, 0., init_nc, k)
        lr = sols[1].u[1]
        push!(lrs, lr)
    end

    p = plot(collect(range(λ_l, λ_h, 100)), [lrs], title = "Leaving Rate vs. λ", labels = "Leaving Rate")
    p
end

# Plot the affinity of a competitive network as a function of the log-volume
function lr_λ_comp(rn, λ_l, λ_h)
    landscape, k = setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, 1., feed_rates_c, comp = true)
    affs = Float64[]
    lrs_A = Float64[]
    lrs_B = Float64[]
    for λ in range(λ_l, λ_h, 100)
        params, sols = volumeOscillation(rn, λ, 0., init_c, k)
        lr_A = sols[1].u[1]
        lr_B = sols[1].u[6]
        aff = A(sols, λ, 0., ΔG_A)
        push!(lrs_A, lr_A)
        push!(lrs_B, lr_B)
        push!(affs, af[1])
    end
    plot(collect(range(λ_l, λ_h, 100)), [lrs_A, lrs_B], title = "Leaving Rates vs. λ", labels = ["A" "B"])
end

function affinity_λ_comp(rn, λ_l, λ_h)
    landscape, k = setRateConstants(rn, ΔG_AF, ΔG_FF, ΔG_A, 1., feed_rates_c, comp = true)
    affs_A = Float64[]
    affs_B = Float64[]
    for λ in range(λ_l, λ_h, 100)
        params, sols = volumeOscillation(rn, λ, 0., init_c, k)
        aff_A = A(sols, λ, 0., ΔG_A)
        aff_B = A(sols, λ, 0., ΔG_A, m_idx = 6)
        push!(affs_A, aff_A[1])
        push!(affs_B, aff_B[1])
    end
    plot(collect(range(λ_l, λ_h, 100)), [affs_A, affs_B], title = "Affinities vs. λ", labels = ["A" "B"])
end

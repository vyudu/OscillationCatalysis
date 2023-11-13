using OscillationCatalysis
println("Testing...")

# Landscapes

@test optimalLandscape(3, 5, 1) == [3, 5, 3, 5, 3, 0]
@test optimalLandscape(4, 3, 4) == [5/4, 3, 5/4, 3, 5/4, 3, 5/4, 0]

# Steady States


# Affinity

@test A(sols, λ0, Δλ, ΔG_A) ≈ 1
@test

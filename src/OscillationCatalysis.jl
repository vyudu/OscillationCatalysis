module OscillationCatalysis

include("CRNOscillation.jl")
export volumeOscillation, setRateConstants, A, A_inst, ratesFromDict 

include("CycleOscillation.jl")
export buildRateMatrix, optimalLandscape, affinity
include("GraphOscillation.jl")
# include("BioCosmos.jl")

include("reactionNetworks.jl")
export oneStep_open, twoStep_open, threeStep_open, autocatalysis_open, competitive_A4B3_open, competitive_A4B2_open, competitive_A4B1_open, closed_competition, oneStep_closed, twoStep_closed, threeStep_closed, autocatalysis_closed, competitive_A4B3_closed, competitive_A4B2_closed, competitive_A4B1_closed

end # module OscillationCatalysis

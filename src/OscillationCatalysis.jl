module OscillationCatalysis

using LinearAlgebra, SymPy, DifferentialEquations, Catalyst, Graphs, MetaGraphs, QuadGK
greet() = print("Hello World!")

include("CRNOscillation.jl")
include("CycleOscillation.jl")
include("GraphOscillation.jl")


end # module OscillationCatalysis

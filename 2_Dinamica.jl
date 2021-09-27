include("src\\Dynamics.jl")

kₘᵢₙ = 0.0
kₘₐₓ = 20*π
N = 10000

G = Grid(kₘᵢₙ, kₘₐₓ, N)

L = Liquid()
ϕ = [0.558]
σ = [1.0]
L.setDistribution(ϕ, σ, phi = true)

S = SF(L, G, HS)

τ, Fs, F, Δζ = dynamics(L, G, S, 7.2, 1e-7, 5, 40)

save_data("DAT\\memoria.dat", [τ Fs F Δζ])


include("src\\StructureFactor.jl")

rho = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

kₘᵢₙ = 0.0
kₘₐₓ = 7*π
N = 1000

G = Grid(kₘᵢₙ, kₘₐₓ, N)

for ρ₀ in rho
	L = Liquid()
	ρ = [ρ₀]
	σ = [1.0]
	L.setDistribution(ρ, σ, phi = false)
	S = SF(L, G, HS; VerletWeis = true, SS = false)
	save_data("DAT\\S"*num2text(ρ₀)*".dat", [G.x S])
end


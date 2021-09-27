include("src\\Asymptotic.jl")

kₘᵢₙ = 0.0
kₘₐₓ = 7*π
N = 1000

G = Grid(kₘᵢₙ, kₘₐₓ, N)
#L = Liquid()
#L.T = 0.01
#L.Soft()
#=
function condicion(ϕ₀)
	println("probando ϕ₀ = ", ϕ₀)
	ϕ = [ϕ₀]
	σ = [1.0]
	L.setDistribution(ϕ, σ, phi = true)
	S = SF(L, G, HS)
	#S = SF(L, G, HS; VerletWeis = true)

	iteraciones, gammas, sistema = Asymptotic(L, G, S, flag = true)
	#return sistema == Fluid()
	return sistema == Glass()
end

ϕ_min = 0.1
ϕ_max = 0.9
δϕ = 0.0001
#ϕᵃ = biseccion(condicion, ϕ_min, ϕ_max, δϕ; flag = false) # Aquiles inicia en ϕ_min, la tortuga está en ϕ_max
ϕᵃ = biseccion(condicion, ϕ_max, ϕ_min, δϕ; flag = false) # Aquiles inicia en ϕ_max, la tortuga está en ϕ_min

print("ϕᵃ = ", ϕᵃ)

=#

Temperaturas = [1e-6]
factor = 2.0

while Temperaturas[end] < 5
	append!(Temperaturas, factor*Temperaturas[end])
end
#println(Temperaturas)

ϕ_min = 0.1
ϕ_max = 1.1
δϕ = 0.0001
	
phi = []
for T in Temperaturas
	println("calculando T =", T)
	L = Liquid()
	L.T = T
	L.Soft()
	function condicion(ϕ₀)
		#println("probando ϕ₀ = ", ϕ₀)
		ϕ = [ϕ₀]
		σ = [1.0]
		L.setDistribution(ϕ, σ, phi = true)
		S = SF(L, G, HS)
		
		iteraciones, gammas, sistema = Asymptotic(L, G, S, flag = false)
		return sistema == Glass()
	end

	ϕᵃ = biseccion(condicion, ϕ_max, ϕ_min, δϕ; flag = false) # Aquiles inicia en ϕ_max, la tortuga está en ϕ_min
	println("ϕᵃ = ", ϕᵃ, "\n")
	append!(phi, ϕᵃ)
end

save_data("DAT\\Arrest_diagram.dat", [phi Temperaturas])


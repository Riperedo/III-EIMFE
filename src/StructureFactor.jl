include("utils.jl")
include("Grid.jl")
include("Liquid.jl")
include("Potentials.jl")

#############################################
# 	Structure Factor for a monodisperese 	#
#	dipole hard spheres colloid				#
#	Wertheim 1971							#
#############################################

#################################
#	Some derivates funtions		#
#################################

# Percus-Yevick solution's coeficients
α₁(ϕ) =  -(1.0 + 2.0*ϕ)^2/(1.0 - ϕ)^4
α₂(ϕ) = 6.0*ϕ*(1.0 + 0.5*ϕ)^2/(1.0 - ϕ)^4
α₃(ϕ) = 0.5*ϕ*α₁(ϕ)

# Some usefull integrals
# ∫(0,1) dx x² j₀(kx)
I₁(k) = k == 0.0 ? 1/3 : (sin(k) - k*cos(k))/k^3
# ∫(0,1) dx x² xj₀(kx)
I₂(k) = k == 0.0 ? 1/4 : (-(k^2 - 2.0)*cos(k) + 2.0*k*sin(k) - 2.0)/k^4
# ∫(0,1) dx x² x³j₀(kx)
I₃(k) = k == 0.0 ? 1/6 : (4.0*k*(k^2 - 6.0)*sin(k) - (k^4 - 12.0*k^2 + 24.0)*cos(k) + 24.0)/k^6

# FT of Wertheim's direct correlation functions
c(ϕ, k) = α₁(ϕ)*I₁(k) + α₂(ϕ)*I₂(k) + α₃(ϕ)*I₃(k)

# Static structure factor
function SF_PY(ϕ :: Real, k :: Real)
	ρ = ϕ2ρ(ϕ)
	return 1.0/(1.0 - 4.0*π*ρ*c(ϕ, k))
end


function SF(L ::  Liquid, G :: Grid, P :: Potential; VerletWeis = false)
	ϕ = L.ϕ[1]
	σ = L.σ[1]
	if VerletWeis & (P == HS)
		σ = (1 - L.ϕ[1]/16)^(1/3)
		ϕ = L.ϕ[1] - (L.ϕ[1]^2)/16
	end
	if L.soft
		T = L.T
		λ = T != 0.0 ? blip(1/T) : 1.0
		λ3 = λ^3
		σ = λ*((1 - λ3*L.ϕ[1]/16)^(1/3))
		ϕ = λ3*L.ϕ[1]*(1 - λ3*L.ϕ[1]/16)
	end
	S = zeros(G.n)
	ρ = L.ρ_total
	for i in 1:G.n
		k = G.x[i]
		Sᵛʷ = SF_PY(ϕ, k*σ)
		U = P.FT(k, L)[1,1]
		Sₛₛ⁻¹ = 1/Sᵛʷ + ρ*U # (1 - ρ(CHS - βU))
		S[i] = 1/Sₛₛ⁻¹
	end
	return S
end


#########################################
#	Structure Factor Blip function 	#
#########################################

function blip(ϵ :: Float64; ν = 6, flag = false)
	λ = 0.0
	dx = 1/1000
	for i in 1:1000
		x = i*dx
		λ -= x*x*exp(-ϵ*(1/x^(2*ν)-2/x^ν+1))
	end
	λ *= 3*dx
	λ += 1.0
	λ = λ^(1/3)
	if flag println("blip function computed ", λ) end
	return λ
end


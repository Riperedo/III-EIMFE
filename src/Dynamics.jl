include("Asymptotic.jl")

#########################
#		Dynamics		#
#########################

function memoria(L, G, S, F, Fs)
	ONE = ones(G.n)
	k = G.x
	integrando = (k.^4).*((S-ONE).^2).*F.*Fs./(S.^2)
	return integral(G, integrando; rule = "trapezoidal")/(6*π*π*L.ρ[1])
end

function decimation(T, Δζ, Fs, F, T_save, Δζ_save, Fs_save, F_save, index)
	N = length(T)
	for n in Int(N//2):N
		if n%2 == 1
			append!(T_save, T[n])
			append!(Δζ_save, Δζ[n])
			append!(Fs_save, Fs[index, n])
			append!(F_save, F[index, n])
		end
	end

	for n in 1:Int(N//2)
		T[n] = T[2*n]
		Δζ[n] = Δζ[2*n]
		Fs[:,n] = Fs[:,2*n]
		F[:,n] = F[:,2*n]
	end
	for n in Int(N//2):N
		T[n] = 0.0
		Δζ[n] = 0.0
		Fs[:,n] .= 0.0
		F[:,n] .= 0.0
	end
	return T, Δζ, Fs, F, T_save, Δζ_save, Fs_save, F_save 
end

function step_free(L, G, S, F, Fs, Δζ, T, dt, n)
	ONE = ones(G.n)
	k = G.x
	α = ONE + dt*(k.*k)./S
	F[:,n+1] = F[:,n]./α
	α = ONE + dt*(k.*k)
	Fs[:,n+1] = Fs[:,n]./α
	Δζ[n+1] = memoria(L, G, S, F[:, n+1], Fs[:, n+1])
	T[n+1] = dt*n
	return T, F, Fs, Δζ
end

function Σᵢ(Δζ, F, n)
	Total = zeros(length(F[:,1]))
	for i in 2:(n-1)
		Total = Total + (Δζ[n+1-i]-Δζ[n-i])*F[:, i]
	end
	return Total
end

function F_dump(L, G, S, F, Δζ, dt, n)
	"""
F_dump(k, n) = α(k)λ(k){Δζₙ₋₁F₁(k) - Σᵢ₌₂ⁿ⁻¹[Δζₙ₊₁₋ᵢ-Δζₙ₋ᵢ]Fᵢ(k)} + Δt⁻¹α(k)Fₙ₋₁(k)
	"""
	ONE = ones(G.n)
	k = G.x
	kc = 2*π*1.305
	λ = ONE./(ONE + (k.*k)/(kc^2))
	α = ONE./(ONE/dt + (k.*k)./S + λ*Δζ[1])
	return α.*(λ.*(Δζ[n-1]*F[:,1] - Σᵢ(Δζ, F, n)) + F[:,n-1]/dt)
end

function F_it(L, G, S, F, Δζ, dt, n)
	"""
F_it(k, n) = α(k)λ(k)Δζₙ[S(k)-F₁(k)]
	"""
	ONE = ones(G.n)
	k = G.x
	kc = 2*π*1.305
	λ = ONE./(ONE + (k.*k)/(kc^2))
	α = ONE./(ONE/dt + (k.*k)./S + λ*Δζ[1])
	return α.*λ.*(S[:]-F[:,1])*Δζ[n]
end

function step(L, G, S, F, Fs, Δζ, T, dt, n)
	"""
Fₙ(k) = α(k)λ(k){Δζₙ₋₁F₁(k) - Σᵢ₌₂ⁿ⁻¹[Δζₙ₊₁₋ᵢ-Δζₙ₋ᵢ]Fᵢ(k)} + Δt⁻¹α(k)Fₙ₋₁(k)
		+ α(k)λ(k)Δζₙ[S(k)-F₁(k)]
	 = F_dump(k, n) + F_it(k, n)
where
F_dump(k, n) = α(k)λ(k){Δζₙ₋₁F₁(k) - Σᵢ₌₂ⁿ⁻¹[Δζₙ₊₁₋ᵢ-Δζₙ₋ᵢ]Fᵢ(k)} + Δt⁻¹α(k)Fₙ₋₁(k)
F_it(k, n) = α(k)λ(k)Δζₙ[S(k)-F₁(k)]
and
α(k) = [Δτ⁻¹I + k²DS⁻¹(k) + λ(k)Δζ₁(k)]⁻¹
	"""
	s = ones(G.n)
	Dump = F_dump(L, G, S, F, Δζ, dt, n)
	dump = F_dump(L, G, s, Fs, Δζ, dt, n)
	Δζ[n+1] = Δζ[n]
	F[:,n+1] = Dump + F_it(L, G, S, F, Δζ, dt, n)
	Fs[:,n+1] = dump + F_it(L, G, s, Fs, Δζ, dt, n)
	Δζ[n+1] = memoria(L, G, S, F[:, n+1], Fs[:, n+1])
	#while true
	for t in 1:10
		F[:,n+1] = Dump + F_it(L, G, S, F, Δζ, dt, n)
		Fs[:,n+1] = dump + F_it(L, G, s, Fs, Δζ, dt, n)
		Δζ_new = memoria(L, G, S, F[:, n+1], Fs[:, n+1])
		if (Δζ[n+1]-Δζ_new)^2 < 1e-10 break end
		Δζ[n+1] = Δζ_new
	end
	T[n+1] = dt*n
	return T, F, Fs, Δζ
end


function dynamics(L, G, S, k_max, dt, nT, decimaciones)
	# grid temporal
	N = 2<<nT
	F = zeros(G.n, N)
	Fs = zeros(G.n, N)
	Δζ = zeros(N)
	T = zeros(N)

	# output
	index = G.find(k_max)
	T_save = []
	Δζ_save = []
	Fs_save = []
	F_save = []
	#initial conditions
	ONE = ones(G.n)
	kc = 2*π*1.305
	F[:, 1] = S
	Fs[:, 1] = ONE
	Δζ[1] = memoria(L, G, S, F[:, 1], Fs[:, 1])
	# first steps
	# free diffusion
	for n in 1:N-1
		T, F, Fs, Δζ = step_free(L, G, S, F, Fs, Δζ, T, dt, n)
	end
	for d in 1:decimaciones
		# decimation
		T, Δζ, Fs, F, T_save, Δζ_save, Fs_save, F_save, = decimation(T, Δζ, Fs, F, T_save, Δζ_save, Fs_save, F_save, index)
		dt *= 2
		for n in Int(N//2):N
			#T, F, Fs, Δζ = step_free(L, G, S, F, Fs, Δζ, T, dt, n-1)
			T, F, Fs, Δζ = step(L, G, S, F, Fs, Δζ, T, dt, n-1)
		end
	end

	# saving final steps
	for n in Int(N//2):N
		append!(T_save, T[n])
		append!(Δζ_save, Δζ[n])
		append!(Fs_save, Fs[index, n])
		append!(F_save, F[index, n])
	end
	
	return T_save, Fs_save, F_save, Δζ_save
end
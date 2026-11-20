# Lennard-Jones energies (full and fractional particles), argon thermal wavelength, and the acceptance rule for λ moves
using Random

function argon_deBroglie(T_σ::Float64)::Float64
    # computes de broglie wavelength for the LJ model of Argon at a given input reduced temperature T_σ = kB T / ϵ

    # -- physical constants (SI)
    h_Js  = 6.62607015e-34      # Planck constant, J s
    kB_J_K = 1.380649e-23        # Boltzmann constant, J/K
    amu_kg = 1.66053906660e-27  # atomic mass unit, kg

    # -- argon LJ parameters (as in Desgranges & Delhommelle 2012)
    ϵ_over_kB_K = 117.05              # epsilon / k_B (K)
    σ_Å = 3.4                       # sigma in Angstroms
    mass_amu = 39.948               # atomic mass of Ar in amu

    # -- convert to SI
    T_K = ϵ_over_kB_K * T_σ             # temperature in K
    σ_m = σ_Å * 1e-10                 # sigma in meters
    m_kg= mass_amu * amu_kg             # mass in kg

    # -- de Broglie wavelength (SI)
    Λ_m = h_Js / sqrt(2*pi * m_kg * kB_J_K * T_K)

    # -- dimensionless quantities in LJ reduced units
    Λ_σ = Λ_m / σ_m


    return (Λ_σ)
end # argon debroglie

function E_12_LJ(rij_squared_σ::Float64)::Float64
    #= Computes the interaction energy between two normal lennard jones particles in LJ units 
    rij_squared_\sigma = squared distance between two particles in lennard jones \sigma =1 units
    =# 
    E_int = (1/rij_squared_σ)^6 - (1/rij_squared_σ)^3
    E_int = 4*E_int 
    return(E_int)
end 

function E_12_frac_LJ(rij_squared_σ::Float64,λ::Int64,λ_max::Int64,ϵ_ξ::Float64,σ_ξ_squared::Float64)::Float64
    #= computes interaction between fractional particle and normal LJ particle according to Eq. 16 of Desgranges & Delhommelle (2012). Note that `M` in Desgranges = λ_max + 1 in our notation
    returns energy in lennard jones units
    # λ and λ_max are not used: the λ dependence enters only through ϵ_ξ and σ_ξ_squared
    # takes ϵ_ξ and σ_ξ_squared precomputed, because they change only when λ changes
    # assumes that rij_squared_σ was computed using the correct minimum image PBC 
    =# 

    E_int = (σ_ξ_squared/rij_squared_σ)^6 - (σ_ξ_squared/rij_squared_σ)^3
    E_int = 4*ϵ_ξ*E_int
    return(E_int)
end #E12_frac_lj 


function potential_1_normal(r_box::Vector{MVector{3,Float64}},ri_box::MVector{3,Float64},i::Int64,r_frac_box::MVector{3,Float64},λ::Int64,λ_max::Int64,N::Int64,L_squared_σ::Float64,r_cut_squared_box::Float64,ϵ_ξ::Float64,σ_ξ_squared::Float64)::Float64
    #= Calculates the sum of the interaction potential between 1 special particle labelled i (e.g. the particle involved in a proposal move) in the r list and all others including the fractional particle
    Inputs:
        -r_box (Vector{MVector{3,Float64}}): positions of the full particles in box units (the first N entries are used)
        -ri_box (static array of Float64s): 3x1 vector of the position of ri you want to use for distances between all other particles in box units
        -i (int): index of particle in r list you want to compute interactions between, used for avoiding self interactions.
        - r_frac_box (static array of Float64s): location of fractional particle in box units
        - λ (int): coupling constant 
        - λ_max (int): largest value of λ, M − 1 in Desgranges & Delhommelle (λ takes the M = λ_max + 1 values 0, …, λ_max)
        - N (int): current number of particles, could be read from _,N = size(r_box) but want to minimize all activity in this loop and we already have to keep track of N outside this function call
        - L_squared_σ (Float64): squared length of box in σ units (LJ this time)
        - r_cut_squared_box  (Float64): cutoff length squared  for the potential in box=1 units
    Returns:
        - E_int_σ (Float64): total interaction energy between particle i at location ri and all other particles including the fractional one in natural units for the potential (i.e. lennard jones units)
    =# 

    E_int_σ = 0

    # first computing interaction with normal particles
    @inbounds for j in 1:N
        if j != i # avoid double counting
            rij_squared_box = euclidean_distance_squared_pbc(ri_box,r_box[j])
            if rij_squared_box < r_cut_squared_box # only evaluate potential if pair is inside cutoff distance
                rij_squared_σ = rij_squared_box  * L_squared_σ   # σ² units
                if rij_squared_σ != 0.0 # r = 0: infinite energy, reject
                    E_int_σ = E_int_σ + E_12_LJ(rij_squared_σ)  
                else
                    return(typemax(Float64))
                end
            end
        end
    end 

    # adding fractional interaction
    rij_squared_box = euclidean_distance_squared_pbc(ri_box,r_frac_box)
    if (rij_squared_box < r_cut_squared_box) && (rij_squared_box != 0.0)# if the distance is zero, E_12_frac_LJ returns NaN
        rij_squared_σ = rij_squared_box  * L_squared_σ  
        E_int_σ = E_int_σ + E_12_frac_LJ(rij_squared_σ,λ,λ_max,ϵ_ξ,σ_ξ_squared)
    elseif rij_squared_box == 0.0 # r = 0: infinite energy, reject
            return(typemax(Float64))
    end
    
    return(E_int_σ)

end # potential_1_normal

function potential_1_frac(r_box::Vector{MVector{3,Float64}},r_frac_box::MVector{3,Float64},λ::Int64,λ_max::Int64,N::Int64,L_squared_σ::Float64,r_cut_squared_box::Float64,ϵ_ξ::Float64,σ_ξ_squared::Float64)::Float64
    #= this function calculates the sum of all interactions between the fractional particle and all the normal particles
        -r_box (Vector{MVector{3,Float64}}): positions of the full particles in box units (the first N entries are used)
        - r_frac_box (static array of Float64s): location of fractional particle in box units
        - λ (int): coupling constant 
        - λ_max (int): largest value of λ, M − 1 in Desgranges & Delhommelle (λ takes the M = λ_max + 1 values 0, …, λ_max)
        - N (int): current number of particles, could be read from _,N = size(r_box) but want to minimize all activity in this loop and we already have to keep track of N outside this function call
        - L_squared_σ (Float64): squared length of box in σ units (LJ this time)
        - r_cut_squared_box  (Float64): cutoff length squared  for the potential in box=1 units
        - ϵ_ξ, σ_ξ_squared (Float64): the fractional particle's LJ parameters at the current λ (precomputed)
    Returns:
        - E_int_σ (Float64): total interaction energy between the fractional particle at r_frac_box and the N full particles, in Lennard-Jones units (typemax(Float64) if the fractional particle sits exactly on a full particle)
    =# 
    E_int_σ = 0

    @inbounds for j in 1:N
        rij_squared_box = euclidean_distance_squared_pbc(r_frac_box,r_box[j])
        if (rij_squared_box < r_cut_squared_box) && (rij_squared_box > 0) # if rij_squared_box ==0 E_12_frac_LJ returns a NaN
            rij_squared_σ = rij_squared_box  * L_squared_σ 
            # no overlap rejection for the fractional particle: at small λ an overlap can still have low energy
            E_int_σ = E_int_σ + E_12_frac_LJ(rij_squared_σ,λ,λ_max,ϵ_ξ,σ_ξ_squared)
        elseif rij_squared_box == 0 # r = 0: infinite energy, reject
            return(typemax(Float64))
        end
    end
    return(E_int_σ)

end # potential_1_frac

function  λ_metropolis_pm1(μ::microstate,
                           μ_prop::microstate, idx_deleted::Int64, # μ_prop: proposed microstate
                           wl::WangLandauVars,sim::SimulationParams)::Bool

        # Acceptance rule for a move Δλ = ±1 (Desgranges & Delhommelle 2012, Eqs. 10–12). Valid only for Δλ = ±1.
        # The form of the criterion and the energies entering ΔE depend on the current and proposed λ and N.
        # assumes that you have already checked that N_proposed is in bounds  (N_min ≤ N_proposed ≤ N_max)

        # first we compute the multiplicative prefactor term involving Q,V,Λ in eqns 10-12-- because the N! terms can cause overflow, we only compute them once we know something about N old vs new
        logQ_diff = wl.logQ_λN[μ.λ+1 , μ.N+1] - wl.logQ_λN[μ_prop.λ+1 , μ_prop.N+1]
        partition_ratio = exp(logQ_diff)
        if (μ.λ > 0 && μ_prop.λ > 0) || (μ.λ == 0 && μ_prop.λ == 0)
            V_Λ_prefactor = sim.V_σ^(μ_prop.N -μ.N) * sim.Λ_σ^(3*μ.N - 3*μ_prop.N)
        elseif μ.λ ==0 && μ_prop.λ > 0
            V_Λ_prefactor = sim.V_σ^(μ_prop.N+1 - μ.N) * sim.Λ_σ^(3*μ.N - 3*(μ_prop.N+1))
        elseif μ.λ > 0 && μ_prop.λ == 0 
           V_Λ_prefactor = sim.V_σ^(μ_prop.N - μ.N - 1) * sim.Λ_σ^(3*(μ.N+1)-3*μ_prop.N)
        end

        # now we compute the exponential part of the criterion having to do with the configurational potential energy
        # and the N-dependent factorial prefactor factorial_prefactor = Nold!/Nnew!

        if μ.N == μ_prop.N # λ changed so change in configurational energy only has to do with fractional particle old vs new
            # Eq. 10 of Desgranges & Delhommelle (2012)
           E_old = potential_1_frac(μ.r_box,μ.r_frac_box,   μ.λ   ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box, μ.ϵ_ξ,μ.σ_ξ_squared ) 
           E_proposed = potential_1_frac(μ_prop.r_box,μ_prop.r_frac_box,      μ_prop.λ    ,sim.λ_max,μ_prop.N,sim.L_squared_σ,sim.r_cut_squared_box, μ_prop.ϵ_ξ,μ_prop.σ_ξ_squared)
           factorial_prefactor = 1
        
        elseif μ.N < μ_prop.N # particle created: λ = λ_max and λ_proposed = 0
            #  change in potential energy comes from fractional particle becoming full particle and the new fractional particle with λ=0 makes no contribution to energy
            # so E_Old = E_old_frac interaction with all others and E_new = E_new_full_particle interaction with all others
            factorial_prefactor = 1/μ_prop.N # only works for ±1
            E_old = potential_1_frac(μ.r_box,μ.r_frac_box,  μ.λ   ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box,  μ.ϵ_ξ,μ.σ_ξ_squared )
            i = μ_prop.N # the new full particle
            E_proposed = potential_1_normal(μ_prop.r_box , μ_prop.r_box[i],i,μ_prop.r_frac_box,μ_prop.λ,sim.λ_max,μ_prop.N,sim.L_squared_σ,sim.r_cut_squared_box, μ_prop.ϵ_ξ,μ_prop.σ_ξ_squared)
                    
        elseif μ.N > μ_prop.N # particle destroyed: λ = 0 and λ_proposed = λ_max
            # old energy is energy of destroyed particle with rest of full particles 
            # new configurational energy is interaction of fractional particle with others
            factorial_prefactor = μ.N # only works for ±1
            E_old = potential_1_normal(μ.r_box,μ.r_box[idx_deleted],idx_deleted,μ.r_frac_box,μ.λ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box,  μ.ϵ_ξ,μ.σ_ξ_squared)
            E_proposed = potential_1_frac(μ_prop.r_box,μ_prop.r_frac_box,  μ_prop.λ   ,sim.λ_max,μ_prop.N ,sim.L_squared_σ,sim.r_cut_squared_box,   μ_prop.ϵ_ξ,μ_prop.σ_ξ_squared)
        end # ΔN logic 

        ΔE = E_proposed - E_old
        exponent = -1*ΔE/sim.T_σ 
        prob_ratio = partition_ratio*V_Λ_prefactor*factorial_prefactor*exp(exponent)
        if prob_ratio > 1
            return(true)
        else
            ζ = rand(sim.rng)
            accept = (prob_ratio > ζ)   #boolean
            return(accept)
        end
end #λ_metropolis_pm1
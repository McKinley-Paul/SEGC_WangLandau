# Structs, initialization, input files and checkpointing (JLD2)
using StaticArrays
using Printf
using JLD2
using Random


mutable struct microstate # holds the variables of the actual current or proposed state
    N::Int64 # number of particles in the microstate
    λ::Int64 # coupling stage of the fractional particle, 0…λ_max
    r_box::Vector{MVector{3,Float64}}  # positions of the full particles (box units), pre-allocated to N_max slots
    r_frac_box::MVector{3,Float64} # position of fractional particle

    # coupling parameters of the fractional particle, recomputed whenever λ changes:
    ϵ_ξ::Float64 # = (λ/M)^(1/3), energy coupling of the fractional particle (LJ units)
    σ_ξ_squared::Float64 # = (σ_ξ)^2 = (λ/M)^(1/2), size coupling of the fractional particle (LJ units)
end

mutable struct SimCache
    # Pre-allocated scratch buffers for the inner loop, to avoid heap allocations. Update them in place (.=).
    ζ_Mvec::MVector{3,Float64} # random displacement
    ζ_idx::Int64 # index of the particle to move
    ri_proposed_box::MVector{3,Float64} # proposed position
    μ_prop::microstate # proposed state for λ moves; mirrors the current microstate between moves
end


struct SimulationParams # fixed for the whole run

    # ~~~~~~~ Input parameters ~~~~~~~~~
    N_max::Int64 # maximum number of atoms in the simulation
    N_min::Int64 # minimum number of atoms in simulation
    T_σ::Float64 # temperature k_B T/ϵ
    Λ_σ::Float64 # thermal de Broglie wavelength in σ units (argon_deBroglie(T_σ) for argon)

    λ_max::Int64 # M − 1 in Desgranges & Delhommelle (2012): λ = 0, …, λ_max, i.e. M stages
    r_cut_σ::Float64 # cutoff to stop evaluating LJ potential

    input_filename::String
    save_directory_path::String # directory for checkpoints, the progress log and the final state
    rng::MersenneTwister # random number generator; pass MersenneTwister(seed) for reproducible runs
    maxiter::Int64 # maximum number of Monte Carlo moves

    # ~~~~~~~ Derived parameters ~~~~~~~~~
    L_σ::Float64 # box length (σ), given directly or read from the input file
    V_σ::Float64 # volume
    L_squared_σ::Float64 # converts squared distances from box units to σ²
    r_cut_box::Float64
    r_cut_squared_box::Float64 # squared cutoff (box units), avoids square roots
    dynamic_δr_max_box::Bool # true: tune δr_max_box during the first stage; false: fixed value passed to init_WangLandauVars

    # inner constructor: computes the derived parameters from the inputs (all keyword arguments)
    function SimulationParams(; N_max,N_min,T_σ,Λ_σ,λ_max,r_cut_σ,input_filename="",L_σ::Float64=0.0,save_directory_path,rng=MersenneTwister(),maxiter=10^9,dynamic_δr_max_box=true)
        # L_σ is given directly or read from the input configuration file
        if input_filename != ""
            _, L_σ, _ = load_configuration(input_filename)
        elseif L_σ == 0.0
            throw(ArgumentError("Must provide either input_filename or L_σ to determine box size"))
        end
        V_σ = L_σ^3
        L_squared_σ = L_σ^2
        r_cut_box = r_cut_σ/L_σ
        r_cut_squared_box = r_cut_box^2

        new(N_max,N_min,T_σ,Λ_σ,λ_max,r_cut_σ,input_filename,save_directory_path,rng,maxiter,
            L_σ,V_σ,L_squared_σ,r_cut_box ,r_cut_squared_box,dynamic_δr_max_box)
    end
end # SimulationParams

mutable struct WangLandauVars
    # Wang–Landau state, updated throughout the run and saved in checkpoints
    logf::Float64 # ln f, the modification factor; starts at 1 (f = e)
    H_λN::Matrix{Int64} # visit histogram, rows λ = 0…λ_max, columns N = 0…N_max; reset at each halving
    logQ_λN::Matrix{Float64} # ln Q(N, V, T, λ) up to a constant, same layout as H_λN
    translation_moves_proposed::Int64
    translation_moves_accepted::Int64
    λ_moves_proposed::Int64
    λ_moves_accepted::Int64
    iters::Int64 # total Monte Carlo moves
    δr_max_box::Float64
    phase2::Bool # unused by the Desgranges schedule; kept for checkpoint compatibility
    flat::Bool # set to true on convergence (ln f < 1e-8)
    
end

function init_WangLandauVars(sim::SimulationParams,δr_max_box::Float64 = typemax(Float64))::WangLandauVars
    if sim.dynamic_δr_max_box == true
        δr_max_box = 0.15/sim.L_σ 
    elseif  ( sim.dynamic_δr_max_box == false ) && (δr_max_box == typemax(Float64) )
        throw("you have sim.dynamic_δr_max_box set to false but did not provide a value of delta r_max_box to init_WangLandauVars")
    end

    logf = 1 
    H_λN = zeros(Int64,sim.λ_max+1,sim.N_max+1)
    logQ_λN=zeros(Float64,sim.λ_max+1,sim.N_max+1) # row λ+1, column N+1
        
    wl = WangLandauVars(logf,H_λN,logQ_λN,0,0,0,0,0,δr_max_box,false,false)
    return(wl)
end



function copy_microstate!(dest::microstate, src::microstate)
    dest.N = src.N
    dest.λ = src.λ
    for i in 1:src.N
        dest.r_box[i] .= src.r_box[i]  # deep copy each MVector's contents
    end
    # copy element-wise: `dest.r_box .= src.r_box` would make the two vectors share the same MVector objects
    dest.r_frac_box .= src.r_frac_box
    dest.ϵ_ξ = src.ϵ_ξ
    dest.σ_ξ_squared = src.σ_ξ_squared
end





function load_configuration(filename::String)::Tuple{Int,Float64,Vector{MVector{3,Float64}}}
    # Reads a configuration file in the format of Allen & Tildesley's initialize.py: line 1 N, line 2 the box length
    # L (σ), then N lines of x y z (σ). Returns N, L_σ and the positions as a vector of MVector{3,Float64} (σ units).
    open(filename, "r") do io
        # Read N (int)
        N = parse(Int, strip(readline(io)))

        # Read L (float)
        L_σ = parse(Float64, strip(readline(io)))

        # Read N lines of 3 floats each
        data = Array{Float64}(undef, N, 3)
        for i in 1:N
            data[i, :] = parse.(Float64, split(readline(io)))
        end

        r_σ = permutedims(data)   # 3×N
        r_σ = matrix_to_svec_positions(r_σ)
        return N, L_σ, r_σ
    end
end # load_configuration(filename)

@inline function matrix_to_svec_positions(r_box::AbstractMatrix{Float64})::Vector{MVector{3,Float64}}
    @assert size(r_box, 1) == 3
    N = size(r_box, 2)

    r = Vector{MVector{3,Float64}}(undef, N)
    @inbounds for j in 1:N
        r[j] = MVector{3,Float64}(
            r_box[1, j],
            r_box[2, j],
            r_box[3, j],
        )
    end
    return r
end

function init_microstate(;sim::SimulationParams,filename::String,λ::Int=0, r_frac_box::MVector{3,Float64}=@MVector [0.0, 0.0, 0.0] )::microstate
    # microstate from a configuration file
    N,L_σ,r_σ = load_configuration(filename)
    r_σ_box = r_σ ./ L_σ # box units, r ∈ [-0.5, 0.5)³
    for i in 1:length(r_σ_box)
        pbc_wrap!(r_σ_box[i])
    end
    # pre-allocate N_max slots so λ-moves can always write to r_box[μ.N] without bounds errors
    r_box = [@MVector zeros(3) for _ in 1:sim.N_max]
    for i in 1:N
        r_box[i] .= r_σ_box[i]
    end

    M = sim.λ_max+1
    ϵ_ξ = (λ/M)^(1/3)
    σ_ξ_squared = (λ/M)^(1/2)

    μstate = microstate(N,λ,r_box,r_frac_box, ϵ_ξ,σ_ξ_squared)
    return(μstate)
end # init microstate

function init_microstate(sim::SimulationParams)::microstate
    # vacuum microstate (N = 0, λ = 0), the anchor Q(N = 0, λ = 0) = 1
    r_box = [@MVector zeros(3) for _ in 1:sim.N_max]
    r_frac_box = @MVector [0.0, 0.0, 0.0]
    return microstate(0, 0, r_box, r_frac_box, 0.0, 0.0)
end # init_microstate from vacuum

@inline function pbc_wrap!(r::MVector{3,Float64})
    @inbounds for i in 1:3
        r[i] = r[i] - round(r[i])
    end
    return nothing
end    



function check_inputs(s::SimulationParams,μ::microstate,wl::WangLandauVars)
    if μ.N > 1
        min_distance,particle_1,particle_2 = min_config_distance(μ.r_box[1:μ.N])
        println("Minimum distance in box units (L_box=1) between particles in the initial configuration is: ",round(min_distance,digits=3), " between particles ",particle_1, " and ", particle_2)
        println("This is compared to length of σ in box units which is: ", round(1/s.L_σ,digits=3))
    else
        println("Starting from N=$(μ.N), skipping minimum distance check.")
    end

    if s.dynamic_δr_max_box == false
        println("SimulationParams.dynamic_δr_max_box is set to false, therefore wl.δr_max_box will not be updated during the run and
        will be set to ", wl.δr_max_box, " for the entire course of the wang landau simulation")
    end
    if (s.dynamic_δr_max_box == false) && (wl.δr_max_box == (0.15/s.L_σ))
        println("wl.δr_max_box: ", wl.δr_max_box)
        throw("You have SimulationParams.dynamic_δr_max_box set to false yet wl.δr_max_box is the
                 default value of 0.15/ L_sigma")
    elseif (s.dynamic_δr_max_box == true) && (wl.δr_max_box != (0.15/s.L_σ))
        println("wl.δr_max_box: ", wl.δr_max_box)
        throw("You have SimulationParams.dynamic_δr_max_box set to true yet wl.δr_max_box is NOT the
                 default value of 0.15/ L_sigma")
    end

    flush(stdout)

end #check_inputs

function init_cache(sim::SimulationParams,initial_μ::microstate)::SimCache
    ζ_Svec = @MVector zeros(3)
    ζ_idx = Int64(0)
    ri_proposed_box = @MVector zeros(3)
    r_box = [@MVector zeros(3) for _ in 1:sim.N_max]
    r_frac_box = @MVector [0.0, 0.0, 0.0]
    μ_prop = microstate(0,
                        0,
                        r_box,
                        r_frac_box,
                        0.0,
                        0.0 )
    copy_microstate!(μ_prop,initial_μ)
    return SimCache(ζ_Svec,ζ_idx,ri_proposed_box,μ_prop)
end


function print_simulation_params(params::SimulationParams,start::Bool=true)
    println()
    if start
        println("           Starting SEGC Wang Landau Simulation        ")
    end
    @printf("T* = %.4f\n", params.T_σ)
    @printf("V = %.4f σ³\n", params.V_σ)
    @printf("Box Length L = %.4f σ\n", params.L_σ)
    @printf("DeBroglie Λ = %.4f σ\n", params.Λ_σ)
    @printf("Max Density (N/V [σ^-3]): ρ = %.4f N\n", params.N_max/params.V_σ)
    
    @printf("Nmax = %d\n", params.N_max)
    @printf("Nmin = %d\n", params.N_min)
    @printf("λmax = %d\n", params.λ_max)

    @printf("r_cut = %.4f σ\n", params.r_cut_σ)
    println("Directory files saved in: ", params.save_directory_path)
end #print_simulation_params

function print_microstate(μ::microstate,print_r::Bool=false)
    println()
    println("Current microstate:")
    @printf("N = %.4f\n", μ.N)
    @printf("λ = %.4f\n", μ.λ)
    if print_r 
        println("Positions of full particles in box units: ")
        display(μ.r_box)

        println("Position of fractional particle:")
        display(μ.r_frac_box)
    end
end #print_microstate

function print_wl(wl::WangLandauVars,verbose=true)
    println()
    println("Wang–Landau variables")
    @printf("logf = %.4f\n", wl.logf)
    @printf("λ_moves_proposed = %f\n", wl.λ_moves_proposed)
    @printf("λ_moves_accepted = %f\n", wl.λ_moves_accepted)
    @printf("λ acceptance ratio = %.4f\n", ( wl.λ_moves_accepted/wl.λ_moves_proposed) )

    @printf("δr_max_box = %.4f\n", wl.δr_max_box)
    @printf("translation_moves_proposed = %f\n", wl.translation_moves_proposed)
    @printf("translation_moves_accepted = %f\n", wl.translation_moves_accepted)
    @printf("translation acceptance ratio = %.4f\n", ( wl.translation_moves_accepted/wl.translation_moves_proposed) )
    @printf("Total monte carlo moves so far = %.1f\n", wl.iters )

    if verbose 
        display("LogQ: ")
        display(wl.logQ_λN)

        display("Histogram λN")
        display(wl.H_λN)
    end
end # print_wl

function initialization_check(sim::SimulationParams, μ::microstate,wl::WangLandauVars)
    check_inputs(sim,μ,wl) # minimum pair distance of the starting configuration, δr_max settings
    print_simulation_params(sim)
    print_microstate(μ,true)
    flush(stdout)
end


function save_wanglandau_jld2(wl::WangLandauVars,sim::SimulationParams, filename::String)
    filepath = joinpath(sim.save_directory_path,filename)
    jldsave(filepath; wl=wl)
end

function load_wanglandau_jld2(filename::String)::WangLandauVars
    return load(filename, "wl")
end

function save_microstate_jld2(μ::microstate,sim::SimulationParams,  filename::String)
    filepath = joinpath(sim.save_directory_path,filename)
    jldsave(filepath; μ=μ)
end

function load_microstate_jld2(filename::String)::microstate
    return load(filename, "μ")
end
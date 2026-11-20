module segc_wl
using StaticArrays
using Dates # timestamps in the progress log
using Statistics

include("initialization.jl")
include("utils.jl")
include("lj.jl")
include("thermo.jl")
# Main loop (run_simulation!), Monte Carlo moves and the Wang–Landau update.

# segc_wl.jl functions:
export run_simulation!, translation_move!,λ_move!,update_wl!, post_run
# initialization functions exports:
export microstate,SimulationParams,init_microstate,check_inputs, print_simulation_params, print_microstate,print_wl, WangLandauVars,init_WangLandauVars, initialization_check, save_wanglandau_jld2, save_microstate_jld2,load_microstate_jld2, load_wanglandau_jld2, load_configuration, SimCache,init_cache,copy_microstate!
# utils exports:
export euclidean_distance, min_config_distance, euclidean_distance_squared_pbc, translate_by_random_vector!, metropolis
# lj exports:
export argon_deBroglie, E_12_LJ, E_12_frac_LJ, potential_1_normal, potential_1_frac, λ_metropolis_pm1
# thermo exports:
export correct_logQ, ideal_gas_logQNVT, ideal_gas_logQ_loggamma


function run_simulation!(sim::SimulationParams, μ::microstate,wl::WangLandauVars,c::SimCache)
    #= Wang–Landau schedule of Desgranges & Delhommelle, J. Chem. Phys. 136, 184107 (2012), Sec. III.C:
         - initial convergence factor f = e (ln f = 1, set in init_WangLandauVars)
         - a stage ends when every (N, l) bin has been visited at least 1000 times (their flatness criterion,
           taken from Ganzenmüller & Camp 2007); then f → √f (ln f halved) and the histogram is reset
         - the run ends when the convergence factor reaches the threshold "f = 10^-8", read as ln f = 10^-8
           (f itself is > 1), as in the original Wang & Landau (2001) criterion
       Moves: 75% translation (of a full or the fractional particle), 25% change in (N, l), with the
       Wang–Landau update after every attempted move. =#
    log_path = joinpath(sim.save_directory_path, "wl_progress_log.txt")
    progress_log = open(log_path,"a") # progress log: one line per halving of ln f
    println(progress_log,"Starting run_simulation!() (Desgranges schedule), time is ",  Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), ", ln f = ", wl.logf)
    flush(progress_log)

    logf_final = 1e-8 # Desgranges 2012: f reduced by square roots "until the minimal threshold value (f = 10^-8)"; read as ln f
    min_visits = 1000 # Desgranges 2012: "all values of N and l have been visited at least 1000 times"

    while (wl.logf ≥ logf_final) && (wl.iters < sim.maxiter)
        ζ = rand(sim.rng)
        if ζ < 0.75                         # propose translational moves 75% of the time
            translation_move!(sim,μ,wl,c)
        else                                # ζ ≥ 0.75: propose (N, l) moves 25% of the time
            λ_move!(sim,μ,wl,c)
        end

        update_wl!(wl,μ)

        if wl.iters % 10_000 == 0 # check the visit criterion every 10^4 moves
            H_min = minimum(wl.H_λN[:,(sim.N_min+1) : (sim.N_max+1) ]) # the zeroth particle sits in wl.H_λN[:,1]
            if H_min ≥ min_visits
                wl.logf = 0.5*wl.logf  # f → √f
                wl.H_λN = zeros(Int64,sim.λ_max+1,sim.N_max+1)
                println(progress_log,"New WL stage, ln f = ",wl.logf,", total iters ", wl.iters, ", time ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
                flush(progress_log)
            end
        end

        if wl.iters % 1_000_000 == 0 # save checkpoint every 10^6 moves
            save_microstate_jld2(μ,sim, "microstate_checkpoint.jld2")
            save_wanglandau_jld2(wl,sim, "wl_checkpoint.jld2")
        end
    end
    wl.logf < logf_final && (wl.flat = true) # converged (the phase2 field is unused; kept so that checkpoint files stay compatible)
    println(progress_log,"Finished: ln f = ", wl.logf, " after ", wl.iters, " total iterations, ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    close(progress_log)
end #run_simulation

function translation_move!(sim::SimulationParams,μ::microstate,wl::WangLandauVars,c::SimCache)
    wl.translation_moves_proposed += 1

    c.ζ_idx = rand(sim.rng,1:(μ.N+1)) # pick a particle to move; index N+1 is the fractional particle
    if c.ζ_idx ≤ μ.N # move normal particle    
        c.ri_proposed_box .= μ.r_box[c.ζ_idx]
        translate_by_random_vector!(c.ri_proposed_box, wl.δr_max_box, sim.rng,c) # trial position (box units)
        pbc_wrap!(c.ri_proposed_box)   # PBC
        E_proposed = potential_1_normal(μ.r_box,c.ri_proposed_box,c.ζ_idx,μ.r_frac_box,μ.λ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box,  μ.ϵ_ξ,μ.σ_ξ_squared ) 

        if E_proposed == typemax(Float64) # check overlap
            nothing # reject the move and recount this state for histogram and partition function 
        else   
            E_old =potential_1_normal(μ.r_box,μ.r_box[c.ζ_idx],c.ζ_idx,μ.r_frac_box,μ.λ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box,  μ.ϵ_ξ,μ.σ_ξ_squared ) 
            ΔE = E_proposed - E_old 
            accept = metropolis(ΔE,sim.T_σ,sim.rng) 
            if accept
                μ.r_box[c.ζ_idx] .= c.ri_proposed_box
                c.μ_prop.r_box[c.ζ_idx] .= c.ri_proposed_box
                wl.translation_moves_accepted += 1
            end
        end

    else # move fractional particle
        c.ri_proposed_box .= μ.r_frac_box
        translate_by_random_vector!(c.ri_proposed_box,wl.δr_max_box,sim.rng,c) # trial position (box units)
        pbc_wrap!(c.ri_proposed_box)   # PBC
        if μ.λ == 0 # always accepted: at λ = 0 the fractional particle does not interact
            μ.r_frac_box .= c.ri_proposed_box
            c.μ_prop.r_frac_box .= c.ri_proposed_box
            wl.translation_moves_accepted += 1
        else 
            E_proposed = potential_1_frac(μ.r_box, c.ri_proposed_box  ,μ.λ ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box,  μ.ϵ_ξ,μ.σ_ξ_squared )

            E_old =potential_1_frac(μ.r_box,μ.r_frac_box,μ.λ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box,  μ.ϵ_ξ,μ.σ_ξ_squared )
            ΔE = E_proposed - E_old 
                accept = metropolis(ΔE,sim.T_σ,sim.rng) 
                if accept
                    μ.r_frac_box .= c.ri_proposed_box
                    c.μ_prop.r_frac_box .= c.ri_proposed_box
                    wl.translation_moves_accepted += 1
                end
        end

    end # full or fractional particle

    if sim.dynamic_δr_max_box == true   
        if wl.logf == 1 # tune δr_max_box during the first Wang–Landau stage only
            if (wl.translation_moves_accepted/wl.translation_moves_proposed > 0.55) && wl.δr_max_box < 1.0 # target ~50% acceptance (Allen & Tildesley, p. 159)
                # the cap δr_max_box < 1 (one box length) stops unbounded growth in dilute systems, where nearly every move is accepted
                wl.δr_max_box = wl.δr_max_box * 1.05

            elseif wl.translation_moves_accepted/wl.translation_moves_proposed  < 0.45
                wl.δr_max_box = wl.δr_max_box*0.95
            end 
        end
    end
end# translation move

function λ_move!(sim::SimulationParams,μ::microstate,wl::WangLandauVars,c::SimCache)
    # Δλ = ±1 moves. Stepping below λ = 0 deletes a particle (N → N−1, λ → λ_max); stepping above λ_max inserts one
    # (N → N+1, λ → 0). Acceptance: λ_metropolis_pm1 (Desgranges & Delhommelle 2012, Eqs. 10–12).
    wl.λ_moves_proposed += 1

    λ_proposed = μ.λ + 2*rand(sim.rng,Bool) - 1 # λ ± 1; −1 or λ_max + 1 signals a change in N
    
    if -1 < λ_proposed ≤ sim.λ_max # only λ changes
        c.μ_prop.λ = λ_proposed
        idx_deleted=0

        c.μ_prop.ϵ_ξ = ( c.μ_prop.λ/(sim.λ_max + 1) )^(1/3)
        c.μ_prop.σ_ξ_squared = (c.μ_prop.λ/(sim.λ_max + 1) )^(1/2) 

        accept = λ_metropolis_pm1(μ, c.μ_prop,idx_deleted, wl,sim)
        if accept == true
            wl.λ_moves_accepted += 1
            μ.λ = c.μ_prop.λ
            μ.ϵ_ξ = c.μ_prop.ϵ_ξ
            μ.σ_ξ_squared = c.μ_prop.σ_ξ_squared
        else # reset the cache if state is rejected 
            c.μ_prop.λ = μ.λ # reset proposed λ to current λ if move rejected
            c.μ_prop.ϵ_ξ = μ.ϵ_ξ 
            c.μ_prop.σ_ξ_squared = μ.σ_ξ_squared
        end

    elseif λ_proposed == -1 # decrement N 
        # a randomly chosen full particle becomes the fractional particle at λ = λ_max
        c.μ_prop.λ = sim.λ_max     
        c.μ_prop.N = μ.N-1
        if c.μ_prop.N < sim.N_min # reject: N would fall below N_min
            c.μ_prop.λ = μ.λ # reset the cache
            c.μ_prop.N += 1
            return(nothing)
        end

        idx_deleted = rand(sim.rng,1:μ.N)
        # swap last particle into deleted slot
        if idx_deleted != μ.N
            c.μ_prop.r_box[idx_deleted] .= μ.r_box[μ.N]
        end
        # deleted particle becomes the new fractional particle: its position must go to r_frac_box.
        # μ.r_box[idx_deleted] is still the deleted position here (the swap above only touched c.μ_prop.r_box)
        c.μ_prop.r_frac_box .= μ.r_box[idx_deleted]

        c.μ_prop.ϵ_ξ = ( c.μ_prop.λ/(sim.λ_max + 1) )^(1/3)
        c.μ_prop.σ_ξ_squared = (c.μ_prop.λ/(sim.λ_max + 1) )^(1/2)

        accept = λ_metropolis_pm1(μ, c.μ_prop,idx_deleted, wl,sim)
        if accept == true
            wl.λ_moves_accepted += 1
            μ.λ = c.μ_prop.λ
            # set r_frac_box to deleted position BEFORE swapping r_box (μ.r_box[idx_deleted] still valid here)
            μ.r_frac_box .= μ.r_box[idx_deleted]
            if idx_deleted != μ.N # only do the swap if we actually deleted a particle that wasn't the last one
                μ.r_box[idx_deleted] .= μ.r_box[μ.N] # swap last particle into deleted slot in actual microstate
            end
            μ.N -= 1

            μ.ϵ_ξ = c.μ_prop.ϵ_ξ
            μ.σ_ξ_squared = c.μ_prop.σ_ξ_squared

        else # reset the cache state if rejected
            c.μ_prop.λ = μ.λ
            c.μ_prop.N += 1
            c.μ_prop.r_box[idx_deleted] .= μ.r_box[idx_deleted] # put back the deleted particle in the cache
            c.μ_prop.r_frac_box .= μ.r_frac_box # restore ghost position in cache
            c.μ_prop.ϵ_ξ = μ.ϵ_ξ
            c.μ_prop.σ_ξ_squared = μ.σ_ξ_squared
        end

    elseif λ_proposed > sim.λ_max # increment N 
        c.μ_prop.λ = 0 
        c.μ_prop.N = μ.N+1
         if c.μ_prop.N  > sim.N_max # reject: N would exceed N_max
            c.μ_prop.λ = μ.λ # reset the cache if move goes out of bounds
            c.μ_prop.N -= 1
            return(nothing)
        end
        c.μ_prop.r_box[c.μ_prop.N] .= μ.r_frac_box # the fractional particle becomes full particle N+1
        @inbounds for i in 1:3
            c.μ_prop.r_frac_box[i] = rand(sim.rng) - 0.5
        end
        idx_deleted = 0

        c.μ_prop.ϵ_ξ = ( c.μ_prop.λ/(sim.λ_max + 1) )^(1/3)
        c.μ_prop.σ_ξ_squared = (c.μ_prop.λ/(sim.λ_max + 1) )^(1/2) 


        accept = λ_metropolis_pm1(μ, c.μ_prop,idx_deleted, wl,sim)
        if accept == true
            wl.λ_moves_accepted += 1
            μ.λ = 0
            μ.N += 1
            μ.r_frac_box .= c.μ_prop.r_frac_box
            μ.r_box[μ.N] .= c.μ_prop.r_box[μ.N] # add the fractional particle to the μ.N position in actual microstate

            μ.ϵ_ξ = c.μ_prop.ϵ_ξ
            μ.σ_ξ_squared = c.μ_prop.σ_ξ_squared 
        else # reset the cache state if rejected
            c.μ_prop.λ = μ.λ
            c.μ_prop.N -= 1
            c.μ_prop.r_frac_box .= μ.r_frac_box
            # c.μ_prop.r_box[μ.N+1] need not be reset: slots above N are never read
            c.μ_prop.ϵ_ξ = μ.ϵ_ξ 
            c.μ_prop.σ_ξ_squared = μ.σ_ξ_squared
        end
    else
        throw("Error in λ move ΔN control flow")
    end # ΔN control flow

end

function update_wl!(wl::WangLandauVars,μ::microstate)
    wl.logQ_λN[μ.λ+1,μ.N+1] += wl.logf
    wl.H_λN[μ.λ+1,μ.N+1] += 1
    wl.iters +=1
end #update_wl!

function post_run(sim::SimulationParams,μ::microstate,wl::WangLandauVars)
    println("Wang Landau converged or reached max iterations, logf has reached ", wl.logf, " and convergence is achieved when logf reaches ", 1e-8)
    println("Iterations: ", (wl.translation_moves_proposed+wl.λ_moves_proposed), " with maxiters: ", sim.maxiter )
    println("Total translation moves proposed: ", wl.translation_moves_proposed, ", translation moves accepted: ", wl.translation_moves_accepted, ", Acceptance ratio: ", wl.translation_moves_accepted/wl.translation_moves_proposed)
    println("Total λ moves proposed: ", wl.λ_moves_proposed, ", λ moves accepted: ", wl.λ_moves_accepted, ", Acceptance ratio: ", wl.λ_moves_accepted/wl.λ_moves_proposed)
    save_wanglandau_jld2(wl,sim,"final_wl.jld2")
    save_microstate_jld2(μ,sim,"final_microstate.jld2")
    println("Final Microstate and Wang Landau variables saved")
end

end
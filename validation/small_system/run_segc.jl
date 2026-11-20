# One converged SEGC-WL run of a small LJ system; writes lnQ.csv (Q(0) = 1) and time.txt to outdir.
# usage (see run_all.sh): julia --project=<SEGC_WangLandau> run_segc.jl seed outdir N_max L_σ T_σ      (r_c = L/2)
using segc_wl, Random, Printf
seed, outdir = parse(Int, ARGS[1]), ARGS[2]; Nmax, L, T = parse(Int, ARGS[3]), parse(Float64, ARGS[4]), parse(Float64, ARGS[5])
mkpath(outdir); cd(outdir); t0 = time()
sim = SimulationParams(N_max = Nmax, N_min = 0, T_σ = T, Λ_σ = argon_deBroglie(T), λ_max = 99, L_σ = L, r_cut_σ = L/2,
                       save_directory_path = outdir, maxiter = 10^14, rng = MersenneTwister(seed))
μ = init_microstate(sim); wl = init_WangLandauVars(sim); c = init_cache(sim, μ)
redirect_stdout(() -> run_simulation!(sim, μ, wl, c), open("stdout.txt", "w"))
lq = correct_logQ(wl); ok = wl.logf < 1e-8
open("lnQ.csv", "w") do io; for (N, v) in enumerate(lq); @printf(io, "%d,%.17g\n", N - 1, v); end; end
write("time.txt", @sprintf("segc seed %d wall_s %.1f iters %d converged %s\n", seed, time() - t0, wl.iters, ok))
foreach(f -> rm(f; force=true), filter(endswith(".jld2"), readdir()))

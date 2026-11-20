# Writes runK/logQ_checkpoint.csv from runK/wl_checkpoint.jld2 for runs that never reached ln f < 1e-8 (the T* = 0.75
# benchmark was stopped before convergence). Same format as logQ_final.csv; the header records ln f
# and the move count at the checkpoint. Run after the jobs are cancelled:
#     julia --project=../.. -t 1 extract_checkpoint.jl            (optional argument: the benchmark directory)
import Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using segc_wl, Printf
dir = length(ARGS) ≥ 1 ? ARGS[1] : @__DIR__
for k in 1:4
    r = joinpath(dir, "run$k"); f = joinpath(r, "wl_checkpoint.jld2")
    isfile(joinpath(r, "DONE")) && (println("run$k  DONE (logQ_final.csv)"); continue)
    isfile(f) || (println("run$k  no checkpoint"); continue)
    wl = try load_wanglandau_jld2(f) catch e; println("run$k  cannot read checkpoint: $e"); continue end
    logQ = correct_logQ(wl)
    open(joinpath(r, "logQ_checkpoint.csv"), "w") do io
        println(io, "# method=SEGC-WL_M100 run=$k T_star=0.75 Lambda_sigma=$(argon_deBroglie(0.75)) V=512 N_max=450 r_cut=3 LRC=not_applied UNCONVERGED logf=$(wl.logf) iters=$(wl.iters)")
        println(io, "N,logQ")
        for N in 0:length(logQ)-1; @printf(io, "%d,%.17g\n", N, logQ[N+1]); end
    end
    @printf("run%d  ln f = %.3e  (%d halvings of 27)  iters = %.3e  min H(N,l) = %d\n", k, wl.logf, round(Int, -log2(wl.logf)), wl.iters, minimum(wl.H_λN))
end

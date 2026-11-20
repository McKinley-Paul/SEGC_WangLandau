# SEGC-WL vs 1/t GC-WL on the same small LJ system: is the mean SEGC ln Q different from the GC-WL mean beyond the
# standard error? Tests both ln Q(N) and the steps d(N) = ln Q(N+1) - ln Q(N) (a biased move shows up in its step).
using Printf, Statistics
dir = ARGS[1]
rd(p) = [parse(Float64, split(l, ",")[2]) for l in readlines(p)]
load(m) = [rd(joinpath(dir, d, "lnQ.csv")) for d in readdir(dir) if startswith(d, m) && isfile(joinpath(dir, d, "lnQ.csv"))]
S, G = load("segc"), load("gcwl")
@printf("%d SEGC runs, %d GC-WL runs\n", length(S), length(G))
function cmp(A, B, label)
    a = reduce(hcat, A); b = reduce(hcat, B)
    ma, mb = vec(mean(a; dims=2)), vec(mean(b; dims=2))
    se = sqrt.(vec(var(a; dims=2)) ./ size(a, 2) .+ vec(var(b; dims=2)) ./ size(b, 2))
    z = (ma .- mb) ./ se
    println("\n", label)
    println("  N    SEGC-GCWL     se      z    sd SEGC  sd GCWL")
    for i in eachindex(z)
        isfinite(z[i]) || continue
        @printf("  %2d   %+8.4f   %7.4f  %+5.1f   %7.4f  %7.4f\n", i - 1, ma[i] - mb[i], se[i], z[i], std(a[i, :]), std(b[i, :]))
    end
    zf = filter(isfinite, z)
    @printf("  max |z| = %.1f over %d values;  count |z| > 2: %d (expect ~%.1f by chance);  count |z| > 3: %d\n",
            maximum(abs.(zf)), length(zf), count(>(2), abs.(zf)), 0.0455 * length(zf), count(>(3), abs.(zf)))
end
cmp(S, G, "ln Q(N) (Q(0) = 1)")
cmp([diff(s) for s in S], [diff(g) for g in G], "step d(N) = ln Q(N+1) - ln Q(N)   (row N)")

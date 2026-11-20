<!-- VS Code preview: ⇧⌘V (Cmd+Shift+V) in this tab, or ⌘K then V for side by side -->

# Benchmark: SEGC-WL against 1/t GC-WL

This benchmark compares the CPU time and accuracy of two methods, each run to its own stopping criterion:

- **SEGC-WL**, this package: M = 100 fractional stages, with the schedule of Desgranges and Delhommelle (2012). ln f
  starts at 1 and is halved once every (N, λ) bin has been visited 1000 times. The run stops at ln f < 10⁻⁸, after
  27 halvings.
- **1/t GC-WL** ([GCWL.jl](https://github.com/McKinley-Paul/GCWL.jl)): whole-particle insertion and deletion with
  the 1/t schedule, stopped by the flatness criterion.

**System.** Lennard-Jones fluid with argon parameters (ε/k_B = 117.05 K, σ = 3.4 Å), V = 512σ³ (L = 8σ),
r_c = 3σ. The standard long-range correction is applied in the analysis. Each run starts from the vacuum, N = 0.
Each method had four independent runs per state point.

**Hardware.** Every run used a single core of the same CPU model, an Intel Xeon 6972P.

**CPU-hours** (`cpu_time.csv`) are wall time × 1 core, summed over all Slurm jobs of a run, because a run continues
from its checkpoint after a time limit. The wall time of each job is the one `main.jl` writes to its timing log at
the end of the job, cross-checked against Slurm accounting (`sacct` Elapsed × AllocCPUS, with AllocCPUS = 1;
`sacct` TotalCPU is reported as 0 on the cluster used). A job stopped by the time limit or cancelled writes no
timing line, so its time is taken from `sacct`. The T* = 0.75 SEGC-WL runs did not converge and were cancelled, so
their CPU-hours are the wall time up to the cancellation: 598.80 h for runs 1 and 2 and 504.30 h for runs 3 and 4.

## Results

ln z_sat is the logarithm of the activity at coexistence, z = Λ⁻³ exp(βμ) in reduced units, as in Desgranges and
Delhommelle (2012). It is found from ln Q(N) by equal weight of the two peaks of P(N). Uncertainties are the sample
standard deviation over the four runs, in the last digits.

| T* | N_max | Method | CPU-hours per run | ln z_sat |
|---|---|---|---|---|
| 1.20 | 350 | 1/t GC-WL | 0.11–0.37 (mean 0.24) | −3.0304(4) |
| 1.20 | 350 | SEGC-WL | 67–93 (mean 78) | −3.0281(14) |
| 0.75 | 450 | 1/t GC-WL | 6.6–8.6 (mean 7.4) | −5.6832(6) |
| 0.75 | 450 | SEGC-WL, **stopped before convergence** | 504–599 (mean 552) | −5.6839(22) |

The TMMC values for the same system are −3.0307 (T* = 1.20) and −5.6831 (T* = 0.75), from the NIST Standard
Reference Simulation Website, LJ SAT-TMMC: liquid–vapor coexistence properties, V = 512σ³, r_c = 3σ with long-range
corrections (https://mmlapps.nist.gov/srs/LJ_PURE/sattmmc.htm).

The T* = 0.75 SEGC-WL runs were stopped before convergence, after 21–24 of the 27 halvings:

| Run | Halvings done | ln f at stop |
|---|---|---|
| 1 | 21 | 4.8×10⁻⁷ |
| 2 | 23 | 1.2×10⁻⁷ |
| 3 | 24 | 6.0×10⁻⁸ |
| 4 | 21 | 4.8×10⁻⁷ |

In these late stages the walk rarely visits the bins at N ≳ 415 in the last fractional stages (λ ≳ 70), and every
bin must reach 1000 visits before ln f is halved. Their ln Q(N) is taken from the last checkpoint
(`extract_checkpoint.jl`).

The SEGC-WL acceptance rules were checked separately against 1/t GC-WL on a small system (V = 64σ³, N_max = 40,
22 and 88 runs; [`../validation/small_system/`](../validation/small_system/README.md)). ln Q(N) and its steps
ln Q(N+1) − ln Q(N) agree to within statistical error at every N.

At T* = 1.20, the mean ln z_sat of the converged SEGC-WL runs deviates from TMMC by +0.088% (3.7σ of the combined
uncertainty), and the mean P_sat by +0.76% (3.1σ).

`results.md` in each folder gives the full output of the analysis: coexistence densities, saturation pressure,
deviations from TMMC and per-run values.

## Layout

```
segc_vs_gcwl_T1.20/              segc_vs_gcwl_T0.75/
  run1–run4/                       run1–run4/
    main.jl                          main.jl                   the run script (resumes from its checkpoint)
    logQ_final.csv                   logQ_checkpoint.csv       ln Q(N), Q(0) = 1, no LRC (T0.75: unconverged)
    wl_progress_log.txt              wl_progress_log.txt       one line per halving of ln f
  cpu_time.csv                     cpu_time.csv              CPU-hours, final ln f and moves per run
  gcwl_reference/                  gcwl_reference/           the 1/t GC-WL runs: logQ_final.csv and cpu_time.csv
  job.sbatch                       job.sbatch                Slurm array job (set account, QOS, CPU constraint)
  results.md                       results.md                analysis output
                                   extract_checkpoint.jl     ln Q from the checkpoints of stopped runs
```

Each `logQ_*.csv` file has a one-line header that records the method, run, T*, Λ and the number of Monte Carlo
moves, followed by `N,logQ` rows.

The analysis script (`analyze.jl`: long-range correction, coexistence, comparison with TMMC) uses functions from
GCWL.jl and is kept in that repository, under `benchmarks/segc_vs_gcwl_T*/`. To regenerate the tables, from the
directory that holds the GCWL.jl clone:

```bash
(cd GCWL.jl/benchmarks/segc_vs_gcwl_T1.20 && julia --project=../.. analyze.jl /path/to/SEGC_WangLandau/benchmarks/segc_vs_gcwl_T1.20)
(cd GCWL.jl/benchmarks/segc_vs_gcwl_T0.75 && julia --project=../.. analyze.jl /path/to/SEGC_WangLandau/benchmarks/segc_vs_gcwl_T0.75)
```

Each writes `results.md` in GCWL.jl's folder, next to `analyze.jl`; the copies here are taken from there.

To rerun a benchmark, edit the `#SBATCH` lines in `job.sbatch` and run `sbatch job.sbatch` from its folder.

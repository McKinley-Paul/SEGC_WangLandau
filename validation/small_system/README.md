<!-- VS Code preview: ⇧⌘V (Cmd+Shift+V) in this tab, or ⌘K then V for side by side -->

# Validation: SEGC-WL against 1/t GC-WL on a small system

This test checks the SEGC-WL acceptance rules against an independent method. Any error in an insertion, deletion,
λ move or crossing between λ stages would bias ln Q(N) in the same direction in every run. Wang–Landau saturation
error, in contrast, changes sign from run to run. The reference method is 1/t grand canonical Wang–Landau
(GCWL.jl), which inserts and deletes whole particles and has no fractional stages.

**System.** Lennard-Jones fluid with the argon thermal wavelength, at T* = 1.2 with V = 64σ³ (L = 4σ). N runs from
0 to 40, so ρ goes up to 0.625σ⁻³. The cutoff is r_c = 2σ, with no long-range correction.

**Runs.** There were 22 SEGC-WL runs (M = 100, ln f halved until ln f < 10⁻⁸; seeds 3001–3022) and 88 1/t GC-WL
runs (seeds 4001–4088). All runs converged. The mean wall time per run was 877 s for SEGC-WL and 6.4 s for GC-WL.

**Test.** For each N, the difference between the two methods' mean ln Q(N) is divided by its standard error. The
same z-test is applied to the steps ln Q(N+1) − ln Q(N), because a biased move appears in its own step.

**Result** (`results.txt`).

- The largest |z| over N = 0–40 is 1.1 for ln Q(N) and 2.0 for the steps. No value has |z| > 2, against about 1.8
  expected by chance.
- The standard error of a step is 0.015–0.035. A bias of about 0.05–0.1 in any single step would therefore have
  been detected; biases ≲ 0.02 per step are not excluded.
- The run-to-run standard deviation of ln Q(N) is 2–9 times larger for SEGC-WL than for GC-WL (0.56 against 0.26
  at N = 40, 0.12 against 0.013 at N = 1).

The two implementations agree to within statistical error at every N.

## Files

| File | Contents |
|---|---|
| `run_segc.jl`, `run_gcwl.jl` | one run of each method |
| `run_all.sh` | all 110 runs and the comparison (about 30 min on 11 cores) |
| `compare.jl` | the z-tests |
| `results.txt` | output of `compare.jl` |
| `runs/*/` | ln Q(N) and timing of each run |

`run_all.sh` needs GCWL.jl. Set `GCWL=/path/to/GCWL.jl`, or clone it next to this repository.

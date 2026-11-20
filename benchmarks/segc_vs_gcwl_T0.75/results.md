<!-- VS Code preview: ⇧⌘V (Cmd+Shift+V) in this tab, or ⌘K then V for side by side -->

# SEGC-WL vs 1/t GC-WL: LJ argon, T* = 0.75, V = 512σ³, N_max = 450, r_c = 3σ + LRC

Mean ± sample standard deviation over the runs of each method. The SEGC-WL runs were stopped before reaching ln f < 1e-8 and use ln Q at their last checkpoint. CPU-hours = wall time × 1 core summed over every job of a run (cpu_time.csv). ln y_sat here is the saturation activity called ln z_sat by NIST and by Desgranges and Delhommelle.

| Method | runs | CPU-hours per run (mean; min–max) | ln y_sat | ρ_liq | ρ_vap | P_sat |
|---|---|---|---|---|---|---|
| 1/t GC-WL | 4 | 7.4; 6.6–8.6 | -5.6832 ± 0.00059 | 0.821709 ± 0.00092 | 0.00363074 ± 3.8e-06 | 0.00263566 ± 2.1e-06 |
| SEGC-WL (M = 100), unconverged | 4 | 551.5; 504.3–598.8 | -5.68389 ± 0.0022 | 0.821858 ± 0.0019 | 0.00360763 ± 3.3e-05 | 0.00261398 ± 4.7e-05 |
| NIST SAT-TMMC | 5 | — | -5.6831 ± 0.0004 | 0.82191 ± 0.0003 | 3.6296e-03 ± 1e-06 | 2.6353e-03 ± 7e-07 |

Deviation of the mean from NIST, relative and in units of the combined uncertainty √(σ_NIST² + sem²):

| Method | ln y_sat | ρ_liq | ρ_vap | P_sat |
|---|---|---|---|---|
| 1/t GC-WL | -0.001% (-0.2σ) | -0.025% (-0.4σ) | +0.032% (0.5σ) | +0.016% (0.3σ) |
| SEGC-WL (M = 100), unconverged | -0.014% (-0.6σ) | -0.006% (-0.1σ) | -0.605% (-1.3σ) | -0.807% (-0.9σ) |

Relative uncertainty (standard deviation over runs / |mean|) and relative error of the mean vs NIST SAT-TMMC, |mean − TMMC|/|TMMC|:

| Method | quantity | ln y_sat | ρ_liq | ρ_vap | P_sat |
|---|---|---|---|---|---|
| 1/t GC-WL | rel. uncertainty | 0.0103% | 0.111% | 0.106% | 0.0795% |
| 1/t GC-WL | rel. error vs TMMC | 0.00144% | 0.0246% | 0.0319% | 0.0157% |
| SEGC-WL (M = 100), unconverged | rel. uncertainty | 0.0394% | 0.232% | 0.918% | 1.8% |
| SEGC-WL (M = 100), unconverged | rel. error vs TMMC | 0.0135% | 0.00646% | 0.605% | 0.807% |
| NIST SAT-TMMC | rel. uncertainty | 0.00758% | 0.0326% | 0.0289% | 0.0251% |

Per run:

| Method | run | CPU-h | ln y_sat | ρ_liq | ρ_vap | P_sat |
|---|---|---|---|---|---|---|
| 1/t GC-WL | 1 | 7.14 | -5.68346 | 0.82296 | 3.63220e-03 | 2.63531e-03 |
| 1/t GC-WL | 2 | 7.11 | -5.68384 | 0.82169 | 3.62499e-03 | 2.63290e-03 |
| 1/t GC-WL | 3 | 8.61 | -5.68247 | 0.82079 | 3.63302e-03 | 2.63674e-03 |
| 1/t GC-WL | 4 | 6.56 | -5.68303 | 0.82139 | 3.63274e-03 | 2.63772e-03 |
| SEGC-WL (M = 100), unconverged | 1 | 598.80 | -5.68499 | 0.81966 | 3.63318e-03 | 2.64235e-03 |
| SEGC-WL (M = 100), unconverged | 2 | 598.80 | -5.68646 | 0.82432 | 3.63177e-03 | 2.65138e-03 |
| SEGC-WL (M = 100), unconverged | 3 | 504.30 | -5.68156 | 0.82156 | 3.60306e-03 | 2.61459e-03 |
| SEGC-WL (M = 100), unconverged | 4 | 504.30 | -5.68254 | 0.82189 | 3.56250e-03 | 2.54759e-03 |

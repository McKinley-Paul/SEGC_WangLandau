<!-- VS Code preview: ⇧⌘V (Cmd+Shift+V) in this tab, or ⌘K then V for side by side -->

# SEGC-WL vs 1/t GC-WL: LJ argon, T* = 1.20, V = 512σ³, N_max = 350, r_c = 3σ + LRC

Mean ± sample standard deviation over the four runs of each method. CPU-hours = wall time × 1 core summed over every job of a run (cpu_time.csv). ln y_sat here is the saturation activity called ln z_sat by NIST and by Desgranges and Delhommelle.

| Method | runs | CPU-hours per run (mean; min–max) | ln y_sat | ρ_liq | ρ_vap | P_sat |
|---|---|---|---|---|---|---|
| 1/t GC-WL | 4 | 0.2; 0.1–0.4 | -3.03039 ± 0.00041 | 0.564661 ± 0.0028 | 0.100447 ± 0.00049 | 0.0773142 ± 0.00011 |
| SEGC-WL (M = 100) | 4 | 78.4; 66.8–93.3 | -3.02806 ± 0.0014 | 0.562482 ± 0.0019 | 0.10237 ± 0.0021 | 0.0778146 ± 0.00037 |
| NIST SAT-TMMC | 5 | — | -3.0307 ± 8e-05 | 0.56316 ± 5e-05 | 1.0034e-01 ± 3e-05 | 7.7225e-02 ± 1e-05 |

Deviation of the mean from NIST, relative and in units of the combined uncertainty √(σ_NIST² + sem²):

| Method | ln y_sat | ρ_liq | ρ_vap | P_sat |
|---|---|---|---|---|
| 1/t GC-WL | +0.011% (1.5σ) | +0.267% (1.1σ) | +0.107% (0.4σ) | +0.115% (1.6σ) |
| SEGC-WL (M = 100) | +0.088% (3.7σ) | -0.120% (-0.7σ) | +2.024% (1.9σ) | +0.763% (3.1σ) |

Relative uncertainty (standard deviation over runs / |mean|) and relative error of the mean vs NIST SAT-TMMC, |mean − TMMC|/|TMMC|:

| Method | quantity | ln y_sat | ρ_liq | ρ_vap | P_sat |
|---|---|---|---|---|---|
| 1/t GC-WL | rel. uncertainty | 0.0136% | 0.502% | 0.486% | 0.139% |
| 1/t GC-WL | rel. error vs TMMC | 0.0109% | 0.267% | 0.107% | 0.115% |
| SEGC-WL (M = 100) | rel. uncertainty | 0.0471% | 0.344% | 2.08% | 0.48% |
| SEGC-WL (M = 100) | rel. error vs TMMC | 0.0878% | 0.12% | 2.02% | 0.763% |
| NIST SAT-TMMC | rel. uncertainty | 0.00251% | 0.00927% | 0.0299% | 0.0184% |

Per run:

| Method | run | CPU-h | ln y_sat | ρ_liq | ρ_vap | P_sat |
|---|---|---|---|---|---|---|
| 1/t GC-WL | 1 | 0.21 | -3.03033 | 0.56471 | 1.00801e-01 | 7.73740e-02 |
| 1/t GC-WL | 2 | 0.37 | -3.03075 | 0.56167 | 9.98141e-02 | 7.72255e-02 |
| 1/t GC-WL | 3 | 0.11 | -3.03065 | 0.56846 | 1.00310e-01 | 7.72214e-02 |
| 1/t GC-WL | 4 | 0.29 | -3.02984 | 0.56381 | 1.00860e-01 | 7.74360e-02 |
| SEGC-WL (M = 100) | 1 | 93.35 | -3.02847 | 0.56386 | 1.04572e-01 | 7.82704e-02 |
| SEGC-WL (M = 100) | 2 | 84.31 | -3.02600 | 0.56138 | 1.02231e-01 | 7.79601e-02 |
| SEGC-WL (M = 100) | 3 | 69.29 | -3.02847 | 0.56435 | 9.95124e-02 | 7.75819e-02 |
| SEGC-WL (M = 100) | 4 | 66.77 | -3.02929 | 0.56033 | 1.03164e-01 | 7.74460e-02 |

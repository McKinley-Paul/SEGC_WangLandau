<!-- VS Code preview: ⇧⌘V (Cmd+Shift+V) in this tab, or ⌘K then V for side by side -->

# SEGC_WangLandau

A Julia implementation of the expanded grand canonical Wang–Landau algorithm (SEGC-WL) of Desgranges and
Delhommelle [1] for the Lennard-Jones fluid. A run gives the canonical partition functions ln Q(N, V, T) for every
N from 0 to N_max at fixed V and T. From these, the grand canonical partition function, and with it coexistence
properties, follow at any chemical potential.

The package was written as the reference implementation of SEGC-WL for a comparison with 1/t grand canonical
Wang–Landau sampling ([GCWL.jl](https://github.com/McKinley-Paul/GCWL.jl)). The benchmark data are in
[`benchmarks/`](benchmarks/README.md).

## Method

The walk moves in (N, λ), where λ = 0, …, λ_max is the coupling stage of one fractional particle, giving
M = λ_max + 1 stages. Ref. [1] uses M = 100. The fractional particle interacts with the full particles through a
Lennard-Jones potential with scaled parameters:

- energy: ε_ξ = (λ/M)^(1/3) ε;
- size: σ_ξ = (λ/M)^(1/4) σ.

A full particle is inserted by growing the fractional particle through all M stages, and deleted by the reverse.

- **Moves.** 75% are translations of a full particle or of the fractional particle, with the Metropolis criterion.
  25% are λ → λ ± 1 moves, with the acceptance rules of Ref. [1], Eqs. 10–12. Stepping above λ_max makes the
  fractional particle a full particle and starts a new one at λ = 0. Stepping below 0 turns a randomly chosen full
  particle into the fractional particle at λ = λ_max.
- **Wang–Landau schedule** (Ref. [1], Sec. III.C).
  - ln Q(N, V, T, λ) is raised by ln f after every attempted move, starting from ln f = 1.
  - ln f is halved, and the histogram reset, once every (N, λ) bin has been visited at least 1000 times. This is
    checked every 10⁴ moves.
  - The run stops when ln f < 10⁻⁸, after 27 halvings.
- **Normalization.** `correct_logQ` returns the λ = 0 row, normalized with Q(N = 0, V, T) = 1. It requires a run
  with N_min = 0; with N_min > 0 the N = 0 bin is never visited and the result is not normalized (no error is
  raised).
- **Box and potential.** Cubic box with periodic boundary conditions, minimum image and a spherical cutoff r_c. The
  standard long-range correction is not applied in the simulation; it can be added to ln Q(N) afterwards.
- **Displacements.** The maximum displacement is tuned towards 50% acceptance during the first stage only, and
  fixed after it.

**Units.** Lengths are in σ, energies in ε, and the temperature is T* = k_BT/ε. Positions are stored in box units,
in [−0.5, 0.5]³. The thermal de Broglie wavelength Λ (in σ) enters the acceptance rules, so ln Q depends on the
mass. `argon_deBroglie(T*)` gives Λ for argon (ε/k_B = 117.05 K, σ = 3.4 Å, m = 39.948 u). The returned ln Q
includes the factor Λ^(−3N). The reduced partition function Q* = Λ^(3N) Q of Ref. [1] is
`logQ[N+1] + 3N*log(Λ_σ)`.

## Installation

Julia ≥ 1.11 is required; the Manifest was resolved with Julia 1.12.1.

```bash
git clone https://github.com/McKinley-Paul/SEGC_WangLandau.git
cd SEGC_WangLandau
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The package module is `segc_wl`.

## Usage

```julia
using segc_wl, Random

T = 1.2
sim = SimulationParams(N_max = 40, N_min = 0, T_σ = T, Λ_σ = argon_deBroglie(T),
                       λ_max = 99, L_σ = 4.0, r_cut_σ = 2.0,
                       save_directory_path = "out", maxiter = 10^13, rng = MersenneTwister(1))
mkpath("out")                         # the output directory must exist
μ     = init_microstate(sim)          # vacuum start, N = 0, λ = 0
wl    = init_WangLandauVars(sim)
cache = init_cache(sim, μ)
run_simulation!(sim, μ, wl, cache)
logQ  = correct_logQ(wl)              # logQ[N+1] = ln Q(N, V, T), Q(0) = 1
```

`examples/lj_argon_small.jl` is this example as a script. It converges in about 15 minutes on one core. A run can
also start from a configuration file (`init_microstate(sim = sim, filename = ...)`), in the format of the
`initialize.py` program of Allen and Tildesley [2]: N on the first line, the box length on the second, then one
`x y z` line per particle, in σ.

**Checkpoints.** `run_simulation!` writes `wl_checkpoint.jld2` and `microstate_checkpoint.jld2` to
`save_directory_path` every 10⁶ moves, and appends one line per halving of ln f to `wl_progress_log.txt`. The
checkpoints hold the Wang–Landau state and the configuration, but not the random number generator, which belongs to
`sim`. To resume, rebuild `sim` with the same parameters and a new seed, so that the continued run does not repeat
the random numbers already used:

```julia
sim = SimulationParams(N_max = 40, N_min = 0, T_σ = T, Λ_σ = argon_deBroglie(T),
                       λ_max = 99, L_σ = 4.0, r_cut_σ = 2.0,
                       save_directory_path = "out", maxiter = 10^13, rng = MersenneTwister(2))   # new seed
wl    = load_wanglandau_jld2(joinpath("out", "wl_checkpoint.jld2"))
μ     = load_microstate_jld2(joinpath("out", "microstate_checkpoint.jld2"))
cache = init_cache(sim, μ)
run_simulation!(sim, μ, wl, cache)
```

The run scripts in `benchmarks/*/run*/main.jl` do this automatically, with a new seed on every restart.

## Tests and validation

```bash
julia --project=. test/runtests.jl
```

`Pkg.test()` (or `] test`) fails on Julia 1.12 for this project with "expected package Test to be registered"; run
the script as above, as CI does.

The tests take about 1 minute. They are of two kinds.

- **Unit tests:**
  - distances under periodic boundaries, random displacements and the Metropolis criterion;
  - the Lennard-Jones and fractional-particle energies, against reference values;
  - acceptance rates of the λ moves, measured against their acceptance rules;
  - the adaptation of the maximum displacement;
  - initialization, and that the cache stays consistent with the microstate through translation, λ and mixed move
    sequences;
  - the detailed-balance identities (prefactor symmetry and energy antisymmetry) of translation, λ = 0 ↔ 1 and
    insertion/deletion moves.
- **Short end-to-end runs.** ln Q(N) is compared with reference values to within 5% of ln Q: at T* = 1 for N = 1
  (exact) and N = 2 (a second-virial estimate integrated to r = ∞), and with the ideal gas at T* = 10⁶ up to
  N = 108. These check the sampler and the normalization rather than the interaction, which contributes only about
  0.02 to ln Q(2) and is not resolved at this tolerance.

The stronger check is the small-system validation in [`validation/small_system/`](validation/small_system/README.md).
It compares SEGC-WL with 1/t GC-WL for a small Lennard-Jones system (V = 64σ³, N_max = 40), using 22 and 88
independent runs. ln Q(N) agrees to within statistical error at every N.

## Benchmark summary

LJ argon, V = 512σ³, r_c = 3σ with long-range correction, four runs per method, one core each of the same CPU model.
Details are in [`benchmarks/README.md`](benchmarks/README.md).

| T* | N_max | Method | CPU-hours per run | ln z_sat |
|---|---|---|---|---|
| 1.20 | 350 | 1/t GC-WL | 0.11–0.37 | −3.0304(4) |
| 1.20 | 350 | SEGC-WL (this package) | 67–93 | −3.0281(14) |
| 0.75 | 450 | 1/t GC-WL | 6.6–8.6 | −5.6832(6) |
| 0.75 | 450 | SEGC-WL (this package), stopped before convergence | 504–599 | −5.6839(22) |

ln z_sat is the logarithm of the activity at coexistence, z = Λ⁻³ exp(βμ) in reduced units, as in Ref. [1]. TMMC
values for the same system are −3.0307 (T* = 1.20) and −5.6831 (T* = 0.75), from NIST [3].

## Repository layout

| Path | Contents |
|---|---|
| `src/` | the package: `segc_wl.jl` (main loop, moves, Wang–Landau update), `lj.jl` (energies, λ-move acceptance), `initialization.jl` (structs, input, checkpoints), `utils.jl`, `thermo.jl` |
| `test/` | unit and regression tests; `mathematica_verification.nb` holds the reference values |
| `examples/` | a small worked example |
| `initial_configs/` | an fcc starting configuration (N = 108, L = 8σ) used by the tests |
| `benchmarks/` | run scripts, ln Q(N), logs and CPU times of the SEGC-WL vs 1/t GC-WL benchmark |
| `validation/` | the small-system comparison with 1/t GC-WL |

## Citation

If you use this code, please cite it (see `CITATION.cff`) together with Ref. [1].

## References

1. C. Desgranges and J. Delhommelle, "Evaluation of the grand-canonical partition function using expanded
   Wang-Landau simulations. I. Thermodynamic properties in the bulk and at the liquid-vapor phase boundary,"
   *J. Chem. Phys.* **136**, 184107 (2012). doi:10.1063/1.4712023
2. M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*, 2nd ed. (Oxford University Press, 2017);
   example programs at https://github.com/Allen-Tildesley/examples.
3. NIST Standard Reference Simulation Website, Lennard-Jones fluid, SAT-TMMC: liquid–vapor coexistence properties
   (V = 512σ³, r_c = 3σ with long-range corrections), https://mmlapps.nist.gov/srs/LJ_PURE/sattmmc.htm.

## License

MIT; see `LICENSE`.

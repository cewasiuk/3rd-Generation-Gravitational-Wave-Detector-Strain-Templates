# Julia waveform port

This folder contains a direct Julia translation of the `chris` Python waveform code.

Files:
- `constants.jl` defines the shared constants and helper functions.
- `annihilation.jl` contains the annihilation strain model.
- `level_transition.jl` contains the level-transition strain model.
- `Project.toml` declares the Julia dependency used by the port.
- `tests_py_jl/` contains scripts and notebooks that generate Python and Julia outputs and compare them.

Notes:
- The translation stays close to the Python structure so outputs are easy to compare.
- Large integer handling differs between Python and Julia, but the current physics routines use small quantum numbers, so the results remain comparable without special-case handling.
- The Lorentzian fit in `constants.jl` is implemented in pure Julia, so no extra fitting package is required.

Quick start:
1. Open a Julia REPL in this folder or run the scripts with `--project=.`.
2. Instantiate the project with `julia --project=. -e 'using Pkg; Pkg.instantiate()'` if you want `SpecialFunctions` available through the project.
3. Include the desired file, for example `include("annihilation.jl")`.

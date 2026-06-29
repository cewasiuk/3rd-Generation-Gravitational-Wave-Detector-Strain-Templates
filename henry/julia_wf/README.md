# Henry Julia Waveform Port

This folder contains a Julia translation of `henry/python_wf/binary_level_transition_functions.py`.
The port lives in `functions.jl` as module `HenryWF`.

The analytic binary-transition waveform path is available through:

```julia
include("functions.jl")
using .HenryWF

h_of_iota = htilde_plus(
    f_grid, M, r, alpha, Omega0, q, m_i, m_f, eta, Gamma_abs;
    use_z_scaling=false,
    numerical_qc=false,
)
h = h_of_iota(0.0)
```

Notes:

- Julia integers can overflow in ordinary `Int` arithmetic where Python silently grows them. The translated `fac` uses `BigInt`, and floating-point callers convert explicitly.
- The Python `find_cf_root` used `iminuit`. This Julia version uses a small derivative-free optimizer so the file can run without extra packages.
- The Python numerical cloud-mass path depends on SciPy `.npz` table loading plus scattered/cubic interpolation. That table-backed `RelScalar`/`MatchedWaveform` model is not enabled in this self-contained Julia file yet. Use `numerical_qc=false` to use the analytical `q_c` path.

## Comparison Harness

The folder `tests_py_jl` contains small scripts that write comparable text files:

```bash
python henry/julia_wf/tests_py_jl/generate_python_test_files.py
julia --compiled-modules=no henry/julia_wf/tests_py_jl/generate_julia_test_files.jl
```

Outputs are written below:

- `henry/julia_wf/tests_py_jl/test_files/python`
- `henry/julia_wf/tests_py_jl/test_files/julia`

Open `compare_test_files.ipynb` to inspect relative differences.

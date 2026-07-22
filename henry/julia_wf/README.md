# Henry Julia Waveform Port

This folder contains a Julia translation of `henry/python_wf/binary_level_transition_functions.py`.
The port lives in `functions.jl` as module `HenryWF`.

The binary-transition waveform path is available through:

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

Use `htilde_cross(...)` with the same arguments to obtain the cross polarization.
Set `numerical_qc=true` to use Henry's table-backed numerical cloud-mass path; this delegates the scalar cloud-mass calculation to the existing Python/SciPy implementation and caches results by parameter tuple.

Notes:

- Julia integers can overflow in ordinary `Int` arithmetic where Python silently grows them. The translated `fac` uses `BigInt`, and floating-point callers convert explicitly.
- The Python `find_cf_root` used `iminuit`. This Julia version uses a small derivative-free optimizer so the file can run without extra packages.
- The numerical cloud-mass path requires the Python environment used by this repository to have NumPy/SciPy available.
- When Julia automatic differentiation reaches the Python-backed cloud-mass lookup, that lookup is evaluated at the primal parameter values and treated as locally constant; the rest of the waveform remains differentiable in Julia.

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

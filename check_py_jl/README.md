# Test Suite for Python/Julia Strain Functions

This directory contains test scripts to verify the consistency between Python and Julia implementations of gravitational wave strain calculations for boson superradiance.

## Files

- `test_python.py` - Python test script that runs calculations and saves outputs
- `test_julia.jl` - Julia test script that runs calculations and compares with Python
- `python.ipynb` - Python notebook for interactive testing
- `julia.ipynb` - Julia notebook for interactive testing

## Usage

### 1. Run Python Tests First

```bash
cd /home/cherrytree/bosonSR_repo_/check_py_jl
python test_python.py
```

This will:
- Run level transition strain calculation
- Run annihilation strain calculation
- Save outputs to `python_level_transition.pkl` and `python_annihilation.pkl`

### 2. Run Julia Tests

```bash
julia test_julia.jl
```

This will:
- Run the same calculations in Julia
- Save outputs to `julia_level_transition.jls` and `julia_annihilation.jls`
- Compare with Python results (if Pickle.jl is installed)
- Print detailed comparison statistics

### 3. Install Pickle.jl (Optional)

To enable automatic comparison of Python and Julia results:

```julia
using Pkg
Pkg.add("Pickle")
```

## What Gets Tested

### Level Transition Strain
- Time evolution of excited and ground state populations
- Strain envelope in time domain
- FFT to frequency domain
- Lorentzian fit to spectral line

**Key outputs compared:**
- Maximum populations (N_e, N_g)
- Maximum strain amplitude
- Transition frequency ω_tr
- Boson mass μ_a
- Lorentzian fit parameters (center frequency f0, width γ)

### Annihilation Strain
- Frequency-domain strain for boson annihilation
- Plus and cross polarizations
- Characteristic strain h_c

**Key outputs compared:**
- Annihilation line frequency
- Gravitational fine-structure constant α
- Maximum strain amplitudes
- Peak frequencies

## Tolerance Levels

- Relative tolerance: 10⁻⁶ (0.0001%)
- Absolute tolerance: 10⁻¹⁰

Small differences may arise from:
- Different FFT implementations
- Different numerical precision in special functions
- Different interpolation methods

## Expected Output

Both tests should produce outputs with relative differences < 10⁻⁶, indicated by ✓ marks in the comparison output.

"""
Test script for Python strain functions.
Runs both level transition and annihilation strain calculations and saves outputs.
"""
import sys
import os
import numpy as np
import json

# Add parent directory to path to import julia_strain_functions
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from julia_strain_functions import iso_gatom_level_tr_strain, iso_gatom_ann_strain

class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder for numpy types."""
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.float64, np.float32, np.float16)):
            return float(obj)
        if isinstance(obj, (np.int64, np.int32, np.int16)):
            return int(obj)
        if isinstance(obj, np.complex128):
            return {"real": obj.real, "imag": obj.imag, "_type": "complex"}
        if isinstance(obj, complex):
            return {"real": obj.real, "imag": obj.imag, "_type": "complex"}
        return super().default(obj)

def save_results(data, filename):
    """Save results to JSON file."""
    json_filename = filename.replace('.pkl', '.json')
    with open(json_filename, 'w') as f:
        json.dump(data, f, cls=NumpyEncoder, indent=2)
    print(f"Saved: {json_filename}")

def test_level_transition():
    """Test level transition strain calculation."""
    print("\n" + "="*60)
    print("Testing Level Transition Strain (Python)")
    print("="*60)
    
    results = iso_gatom_level_tr_strain(
        M_solar=1e-6,
        a_spin=0.999999,
        alpha=1.0,
        ne=6,
        ng=5,
        m=None,
        distance_kpc=10.0,
        N_e0=1.0,
        N_g0=1.0,
        n_time=10000,  # Reduced for faster testing
        n_fft=2**18,   # Reduced for faster testing
        n_top=400,
        verbose=True
    )
    
    # Extract key results for comparison
    output = {
        'time_yr': results['time_yr'],
        'h_t': results['h_t'],
        'N_e': results['N_e'],
        'N_g': results['N_g'],
        'f_pos': results['f_pos'],
        'H_pos': results['H_pos'],
        'L_full': results['L_full'],
        'lorentz_params': results['lorentz_params'],
        'omega_tr_GeV': results['omega_tr_GeV'],
        'Mu_a': results['Mu_a'],
        # Summary statistics for easier comparison
        'max_Ne': np.max(results['N_e']),
        'max_Ng': np.max(results['N_g']),
        'max_h': np.max(np.abs(results['h_t'])),
        'max_freq_pos': results['f_pos'][np.argmax(np.abs(results['H_pos']))],
        'A': results['lorentz_params'][0],
        'f0': results['lorentz_params'][1],
        'gamma': results['lorentz_params'][2],
        'C': results['lorentz_params'][3],
    }
    
    save_results(output, 'python_level_transition.pkl')
    
    print(f"\nSummary Statistics:")
    print(f"  Max N_e: {output['max_Ne']:.6e}")
    print(f"  Max N_g: {output['max_Ng']:.6e}")
    print(f"  Max h(t): {output['max_h']:.6e}")
    print(f"  Lorentzian center f0: {output['f0']:.6e} Hz")
    print(f"  Lorentzian width γ: {output['gamma']:.6e} Hz")
    
    return output

def test_annihilation():
    """Test annihilation strain calculation."""
    print("\n" + "="*60)
    print("Testing Annihilation Strain (Python)")
    print("="*60)
    
    results = iso_gatom_ann_strain(
        M_solar=3.1e-4,
        mua=2e-16,
        n=4,
        l=None,
        alpha=None,
        distance_kpc=1.0,
        iota=0.0,
        phase=0.0,
        f_min_Hz=1e9,
        f_max_Hz=1e11,
        n_f=5000,  # Reduced for faster testing
        verbose=True
    )
    
    # Extract key results for comparison
    output = {
        'f_Hz': results['f_Hz'],
        'h_plus': results['h_plus'],
        'h_cross': results['h_cross'],
        'h_c': results['h_c'],
        'f_line_Hz': results['f_line_Hz'],
        'alpha': results['alpha'],
        'mua_GeV': results['mua_GeV'],
        'M_solar': results['M_solar'],
        'distance_kpc': results['distance_kpc'],
        'n': results['n'],
        'l': results['l'],
        'iota': results['iota'],
        'phase': results['phase'],
        # Summary statistics
        'max_h_plus': np.max(np.abs(results['h_plus'])),
        'max_h_cross': np.max(np.abs(results['h_cross'])),
        'max_h_c': np.max(results['h_c']),
        'peak_freq': results['f_Hz'][np.argmax(results['h_c'])],
    }
    
    save_results(output, 'python_annihilation.pkl')
    
    print(f"\nSummary Statistics:")
    print(f"  Line frequency: {output['f_line_Hz']:.6e} Hz")
    print(f"  Max |h_plus|: {output['max_h_plus']:.6e}")
    print(f"  Max |h_cross|: {output['max_h_cross']:.6e}")
    print(f"  Max h_c: {output['max_h_c']:.6e}")
    print(f"  Peak frequency: {output['peak_freq']:.6e} Hz")
    
    return output

if __name__ == "__main__":
    print("Starting Python test suite...")
    
    # Test level transition
    level_results = test_level_transition()
    
    # Test annihilation
    ann_results = test_annihilation()
    
    print("\n" + "="*60)
    print("Python tests completed successfully!")
    print("="*60)

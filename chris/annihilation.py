"""
annihilation.py
---------------
Analytic frequency-domain strain for boson-cloud annihilation lines,
plus the private helpers specific to this calculation.

Depends on: constants.py
"""

import numpy as np
import mpmath as mp

from constants import (
    G, M_sun, kpc_to_GeV,
    Hz_to_GeV, GeV_to_Hz, h_GeV_s,
    r_g,
)


# ============================================================
# PRIVATE HELPERS  (annihilation-specific)
# ============================================================

def _mp_E1_vec(z):
    """Vectorised exponential integral E1(z) via mpmath."""
    z = np.asarray(z, dtype=complex)
    return np.vectorize(lambda zz: complex(mp.e1(zz)))(z)


def omega_ann(mua, alpha, n):
    """
    Annihilation line angular frequency  [GeV].

    ω_ann ≈ 2 μ_a (1 - α^2 / (2 n^2))

    Parameters
    ----------
    mua   : float — boson mass μ_a [GeV]
    alpha : float — gravitational fine-structure α = G M μ
    n     : int   — principal quantum number
    """
    return 2.0 * mua * (1.0 - alpha**2 / (2.0 * n**2))


def gamma_ann(l, alpha, M_GeV):
    """
    Annihilation decay rate Γ_a  [GeV].

    Parameters
    ----------
    l     : int   — orbital angular momentum quantum number
    alpha : float — gravitational fine-structure α = G M μ
    M_GeV : float — BH mass [GeV]
    """
    p = 17 if l == 1 else 4 * l + 1
    return G * 1e-10 / r_g(M_GeV)**3 * (
        ((alpha / l) * 0.5)**p + ((alpha / l) * 0.5)**(p + 1)
    )


def h_ann(delta_f_Hz, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV):
    """
    Core annihilation line shape h̃(Δf) — one complex amplitude  [strain/Hz].

    Parameters
    ----------
    delta_f_Hz  : array_like — frequency offset from line centre [Hz]
    mua         : float      — boson mass μ_a [GeV]
    M_solar     : float      — BH mass [solar masses]
    n, l        : int        — principal and orbital quantum numbers
    alpha       : float      — α = G M μ (dimensionless)
    r_kpc       : float      — source distance [kpc]
    omega_a_GeV : float      — annihilation line energy [GeV]
    """
    M_GeV  = M_solar * M_sun
    r      = r_kpc * kpc_to_GeV              # GeV^-1
    Gamma  = gamma_ann(l, alpha, M_GeV)      # GeV
    N_max  = 10.0**76 * (M_solar / 10.0)**2

    delta_E = np.asarray(delta_f_Hz, dtype=float) * Hz_to_GeV
    z       = 1j * delta_E / (Gamma * N_max)

    pref = (1.0 / (2.0 * np.pi)) * np.sqrt(
        4.0 * G / (Gamma * r**2 * omega_a_GeV)
    )

    Ei_part = _mp_E1_vec(z)
    Ei      = np.exp(1j * delta_E / (Gamma * N_max)) * Ei_part

    return pref * np.exp(z) * Ei * h_GeV_s


def h_pcr_ann(f_grid_Hz, mua, M_solar, n, l, alpha, r_kpc, iota, phase):
    """
    Frequency-domain annihilation strain: plus and cross polarisations.

    Parameters
    ----------
    f_grid_Hz : array_like — detector frequency grid [Hz]
    mua       : float      — boson mass μ_a [GeV]
    M_solar   : float      — BH mass [solar masses]
    n, l      : int        — principal and orbital quantum numbers
    alpha     : float      — α = G M μ (dimensionless)
    r_kpc     : float      — source distance [kpc]
    iota      : float      — inclination angle [rad]
    phase     : float      — overall phase offset [rad]

    Returns
    -------
    h_plus    : ndarray — h̃_+(f) [strain/Hz]
    h_cross   : ndarray — h̃_×(f) [strain/Hz]
    f_grid_Hz : ndarray — input frequency grid [Hz]
    f_line_Hz : float   — annihilation line frequency [Hz]
    """
    omega_a_GeV = omega_ann(mua, alpha, n)
    f_a_Hz      = omega_a_GeV * GeV_to_Hz

    f_grid_Hz     = np.asarray(f_grid_Hz, dtype=float)
    delta_f_minus = f_grid_Hz - f_a_Hz
    delta_f_plus  = -(f_grid_Hz - f_a_Hz)

    hmin = h_ann(delta_f_minus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)
    hp   = h_ann(delta_f_plus,  mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)

    c       = np.cos(iota)
    h_plus  = (1.0 + c**2) / 4.0 * (
        np.exp( 1j * phase) * hmin + np.exp(-1j * phase) * hp
    )
    h_cross = (c / (2.0j)) * (
        np.exp( 1j * phase) * hmin - np.exp(-1j * phase) * hp
    )

    return h_plus, h_cross, f_grid_Hz, f_a_Hz


# ============================================================
# MASTER FUNCTION
# ============================================================

def iso_gatom_ann_strain(
    M_solar      = 3.1e-4,
    mua          = 2e-16,
    n            = 4,
    l            = None,
    alpha        = None,
    distance_kpc = 1.0,
    iota         = 0.0,
    phase        = 0.0,
    f_min_Hz     = 1.0,
    f_max_Hz     = 1.0e12,
    n_f          = 50000,
    verbose      = True,
):
    """
    Analytic frequency-domain strain for annihilations (no time-domain
    evolution required).

    Parameters
    ----------
    M_solar      : float      — BH mass [solar masses]
    mua          : float      — boson mass μ_a [GeV]
    n            : int        — principal quantum number (e.g. 4 for 4ℓ)
    l            : int or None — orbital quantum number; default l = n - 1
    alpha        : float or None — α = G M μ; computed from M_solar & mua if None
    distance_kpc : float      — source distance [kpc]
    iota         : float      — inclination angle [rad]
    phase        : float      — overall phase offset [rad]
    f_min_Hz     : float      — minimum frequency [Hz]
    f_max_Hz     : float      — maximum frequency [Hz]
    n_f          : int        — number of log-spaced frequency samples
    verbose      : bool       — print diagnostics if True

    Returns
    -------
    dict with keys:
        'f_Hz'         : frequency grid [Hz]
        'h_plus'       : h̃_+(f) [strain/Hz]
        'h_cross'      : h̃_×(f) [strain/Hz]
        'h_c'          : characteristic strain h_c(f)
        'f_line_Hz'    : annihilation line frequency [Hz]
        'alpha'        : gravitational fine-structure α
        'mua_GeV'      : μ_a [GeV]
        'M_solar'      : BH mass [solar masses]
        'distance_kpc' : source distance [kpc]
        'n', 'l'       : quantum numbers
        'iota'         : inclination [rad]
        'phase'        : phase offset [rad]
    """
    if l is None:
        l = n - 1

    if alpha is None:
        alpha = G * (M_solar * M_sun) * mua

    f_grid_Hz = np.logspace(np.log10(f_min_Hz), np.log10(f_max_Hz), n_f)

    h_plus, h_cross, f_grid_Hz, f_line_Hz = h_pcr_ann(
        f_grid_Hz, mua, M_solar, n, l, alpha, distance_kpc, iota, phase,
    )

    h_c = 2.0 * f_grid_Hz * np.sqrt(np.abs(h_plus)**2 + np.abs(h_cross)**2)

    if verbose:
        print(f"alpha           ≈ {alpha:.3e}")
        print(f"r_g(M)          ≈ {r_g(M_solar * M_sun):.4e} GeV^-1")
        print(f"Annihilation f_ann ≈ {f_line_Hz:.3e} Hz")
        print(f"max |h_plus|    ≈ {np.max(np.abs(h_plus)):.3e}")
        print(f"max h_c         ≈ {np.max(h_c):.3e}")

    return {
        "f_Hz":         f_grid_Hz,
        "h_plus":       h_plus,
        "h_cross":      h_cross,
        "h_c":          h_c,
        "f_line_Hz":    f_line_Hz,
        "alpha":        alpha,
        "mua_GeV":      mua,
        "M_solar":      M_solar,
        "distance_kpc": distance_kpc,
        "n":            n,
        "l":            l,
        "iota":         iota,
        "phase":        phase,
    }

"""
level_transition.py
-------------------
Master function for the gravitational-wave strain from a
boson-cloud level transition (e.g. 6g → 5g).

Depends on: constants.py
"""

import numpy as np

from constants import (
    G, M_sun, kpc_to_GeV, GeV_to_yrinv,
    r_g, r_plus, Omega_plus,
    super_gamma, gamma_t_6g_to_5g,
    dNedt, dNgdt, strain_envelope,
    fit_lorentzian_to_fft,
)


def iso_gatom_level_tr_strain(
    M_solar      = 1e-6,
    a_spin       = 0.999999,
    alpha        = 1.0,
    ne           = 6,
    ng           = 5,
    m            = None,
    distance_kpc = 10.0,
    N_e0         = 1.0,
    N_g0         = 1.0,
    n_time       = 100000,
    n_fft        = 2**20,
    n_top        = 400,
    verbose      = True,
):
    """
    Evolve the level populations and strain envelope, then compute
    the frequency-domain strain and fit a Lorentzian to its magnitude.

    Parameters
    ----------
    M_solar      : float — BH mass [solar masses]
    a_spin       : float — dimensionless BH spin (0 ≤ a ≤ 1)
    alpha        : float — gravitational fine-structure constant α = G M μ
    ne           : int   — excited-state principal quantum number (e.g. 6)
    ng           : int   — ground-state principal quantum number (e.g. 5)
    m            : int or None — azimuthal quantum number; default m = ℓ = ng-1
    distance_kpc : float — source distance [kpc]
    N_e0         : float — initial excited-state population
    N_g0         : float — initial ground-state population
    n_time       : int   — number of log-spaced time steps
    n_fft        : int   — FFT length
    n_top        : int   — FFT points used for Lorentzian fit
    verbose      : bool  — print diagnostics if True

    Returns
    -------
    dict with keys:
        'time_yr'        : time grid [yr]
        'h_t'            : strain envelope h(t) (dimensionless)
        'N_e'            : excited-state population history
        'N_g'            : ground-state population history
        'f_pos'          : positive frequency grid [Hz]
        'H_pos'          : FFT of h(t) on f_pos (complex)
        'L_full'         : best-fit Lorentzian (no floor) on f_pos
        'lorentz_params' : (A, f0, gamma, C) from the fit
        'omega_tr_GeV'   : transition angular frequency [GeV]
        'Mu_a'           : boson mass [GeV]
    """
    # Quantum numbers
    l = ng - 1
    if m is None:
        m = l

    # BH mass and boson mass
    M    = M_solar * M_sun
    mu_a = alpha / (G * M)

    # Superradiance rates [GeV] → [1/yr]
    gamma_sre_GeV = super_gamma(ne, l, m, mu_a, M, a_spin)
    gamma_srg_GeV = super_gamma(ng, l, m, mu_a, M, a_spin)
    gamma_sre_yr  = gamma_sre_GeV * GeV_to_yrinv
    gamma_srg_yr  = gamma_srg_GeV * GeV_to_yrinv

    # Transition angular frequency and rate
    omega_tr       = 0.5 * mu_a * alpha**2 * ((1.0 / ng**2) - (1.0 / ne**2))
    gamma_t_ne_GeV = gamma_t_6g_to_5g(alpha, omega_tr, M)
    gamma_t_ne_yr  = gamma_t_ne_GeV * GeV_to_yrinv

    print(f"\n alpha    = {alpha}")
    print(f" omega_tr = {omega_tr}")
    print(f" r_g(M)   = {r_g(M):.4e} GeV^-1")
    print(f" M        = {M:.4e} GeV")

    # Distance in natural units
    dist_GeV_inv = distance_kpc * kpc_to_GeV

    # Time stepping (log-spaced years)
    dt_arr   = np.logspace(-8, 1, num=n_time)
    yr       = 0.0
    Ne, Ng   = N_e0, N_g0
    omega_sr = mu_a

    time_yr  = []
    N_e_hist = []
    N_g_hist = []
    h_hist   = []

    for j in dt_arr:
        r_h     = r_plus(M, a_spin)
        Omega_H = Omega_plus(a_spin, r_h)

        if m * Omega_H > omega_sr:
            dNe = dNedt(Ne, Ng, gamma_sre_yr, gamma_t_ne_yr)
            dNg = dNgdt(Ne, Ng, gamma_srg_yr, gamma_t_ne_yr)

            Ne_next = Ne + dNe * j
            Ng_next = Ng + dNg * j

            if Ne_next < 0 or Ng_next < 0:
                break

            Ne, Ng = Ne_next, Ng_next
            yr    += j

            h_val = strain_envelope(
                dist_GeV_inv, alpha, ne, ng,
                mu_a, gamma_t_ne_GeV, Ng, Ne,
            )

            N_e_hist.append(Ne)
            N_g_hist.append(Ng)
            h_hist.append(h_val)
            time_yr.append(yr)
        else:
            break

    time_yr  = np.array(time_yr)
    N_e_hist = np.array(N_e_hist)
    N_g_hist = np.array(N_g_hist)
    h_hist   = np.array(h_hist, dtype=float)

    if verbose:
        print(f"Max N_e      ≈ {N_e_hist.max():.3e}")
        print(f"Max N_g      ≈ {N_g_hist.max():.3e}")
        print(f"Max h(t)     ≈ {np.max(np.abs(h_hist)):.3e}")
        print(f"omega_tr     ≈ {omega_tr:.3e} GeV")
        print(f"gamma_sre_yr ≈ {gamma_sre_yr:.3e} 1/yr")
        print(f"gamma_srg_yr ≈ {gamma_srg_yr:.3e} 1/yr")
        print(f"gamma_t_yr   ≈ {gamma_t_ne_yr:.3e} 1/yr")
        print(f"gamma_srg    ≈ {gamma_srg_GeV:.3e} 1/GeV")
        print(f"Mu_a         ≈ {mu_a:.3e} GeV")

    f_pos, H_pos, L_full, lorentz_params = fit_lorentzian_to_fft(
        time_yr, h_hist, n_fft=n_fft, n_top=n_top,
    )

    return {
        "time_yr":        time_yr,
        "h_t":            h_hist,
        "N_e":            N_e_hist,
        "N_g":            N_g_hist,
        "f_pos":          f_pos,
        "H_pos":          H_pos,
        "L_full":         L_full,
        "lorentz_params": lorentz_params,
        "omega_tr_GeV":   omega_tr,
        "Mu_a":           mu_a,
    }

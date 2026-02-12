import numpy as np
import math as ma
import numpy as np
import mpmath as mp
from scipy.optimize import curve_fit

# ============================================================
# GLOBAL CONSTANTS
# ============================================================
# Mp        : Planck mass [GeV]
# G         : Newton's constant [GeV^-2]
# hbar_GeV_s: reduced Planck constant ħ [GeV·s]
# s_per_yr  : seconds per year [s]
# GeV_to_yrinv: conversion factor from a rate in GeV to a rate in 1/yr
#
# M_sun     : solar mass [GeV]
# Mpc_to_GeV: 1 Mpc in natural units [GeV]
# kpc_to_GeV: 1 kpc in natural units [GeV]
# S_PER_YR  : seconds per year [s]
Mp  = 1.220890e19                 # GeV
G   = 1.0 / Mp**2                 # GeV^-2
hbar_GeV_s = 6.582119569e-25      # GeV·s
s_per_yr   = 3.154e7
GeV_to_yrinv = (1.0 / hbar_GeV_s) * s_per_yr

M_sun     = 1.98847e30 * 5.60958885e26  # GeV
Mpc_to_GeV = 1.56e38
kpc_to_GeV = Mpc_to_GeV / 1000.0
S_PER_YR   = s_per_yr

h_GeV_s    = 4.1357e-24           # h in GeV·s (E = h f)
Hz_to_GeV  = h_GeV_s   
GeV_to_Hz  = 1.0 / h_GeV_s


# ============================================================
# Isolated Level Transition helpers (frequency-domain analytic strain)
# ============================================================
def r_plus(M, a):
    """Horizon radius r_+ = G M (1 + sqrt(1 - a^2))."""
    return G * M * (1.0 + np.sqrt(1.0 - a**2))   # GeV^-1


def Omega_plus(a, r):
    """Horizon angular frequency Ω_+ = a / (2 r_+)."""
    return a / (2.0 * r)                         # GeV


def omega_bound(mu, alpha, n):
    """Hydrogenic bound-state frequency approximation."""
    return mu * (1.0 - alpha**2 / (2.0 * n**2))


def glm(a, m, l, r, omega):
    """g_{lm} product factor for the small-alpha rate."""
    factor = (a * m - 2.0 * r * omega)**2
    g = 1.0
    for k in range(1, l + 1):
        g *= (k**2 * (1.0 - a**2) + factor)
    return g


# ============================================================
# superradiance rates
# ============================================================
def super_gamma_lowalpha(n, l, m, mu, M, a):
    """Small-alpha superradiance rate Γ_low in GeV."""
    alpha = G * M * mu
    r_p   = r_plus(M, a)
    Omega_p = Omega_plus(a, r_p)
    omega_n = omega_bound(mu, alpha, n)

    C = (
        2.0**(4*l + 1)
        * ma.factorial(n + l)
        / (n**(2*l + 4) * ma.factorial(n - l - 1))
        * (ma.factorial(l) / (ma.factorial(2*l) * ma.factorial(2*l + 1)))**2
    )
    g_lm = glm(a, m, l, r_p, omega_n)

    gamma = (
        2.0 * r_p / M
        * C
        * g_lm
        * (m * Omega_p - omega_n)
        * (alpha**(4*l + 5))
        / G / 4.0
    )
    return gamma


def gamma_wkb(mu, M, a, kappa=3.7):
    """WKB superradiance rate Γ_WKB in GeV."""
    alpha = G * M * mu
    r_g   = G * M
    pref  = 1.0e-7 / r_g
    return pref * np.exp(-kappa * alpha)


def super_gamma(n, l, m, mu, M, a,
                alpha_low=0.1, alpha_high=15.0):
    """Matched superradiance rate (low-alpha + WKB) in GeV."""
    alpha = G * M * mu
    kappa = (l + 0.5)

    if alpha <= alpha_low:
        return super_gamma_lowalpha(n, l, m, mu, M, a)
    if alpha >= alpha_high:
        return gamma_wkb(mu, M, a, kappa=kappa)

    gamma_low = super_gamma_lowalpha(n, l, m, mu, M, a)
    gamma_hig = gamma_wkb(mu, M, a, kappa=kappa)
    log_low   = np.log(gamma_low)
    log_hig   = np.log(gamma_hig)
    w = (alpha - alpha_low) / (alpha_high - alpha_low)
    return np.exp((1.0 - w)*log_low + w*log_hig)


# ============================================================
# transition rate 6g -> 5g
# ============================================================
def gamma_t_6g_to_5g(alpha, omega_tr, M):
    """
    Transition rate Γ_t for 6g → 5g (ℓ = m = 4) in GeV.
    alpha: dimensionless coupling G M mu
    omega_tr: transition angular frequency in GeV
    """
    r_g = G * M
    C = (2**28 * 3**4 * 5**5) / (11**22 * np.pi)
    C *= (32*np.pi/15.0)   # ∫ sin^4θ dΩ
    P_t = C * (G * alpha**12) / (r_g**4)
    return P_t / omega_tr


# ============================================================
# population ODEs and strain envelope
# ============================================================
def dNedt(Ne, Ng, gamma_sr_yr, gamma_t_yr):
    """
    Time derivative dN_e/dt in 1/yr units.

    Parameters
    ----------
    Ne : float
        Excited-state population N_e (dimensionless count).
    Ng : float
        Ground-state population N_g (dimensionless count).
    gamma_sr_yr : float
        Superradiance growth rate for the excited state [1/yr].
    gamma_t_yr : float
        Transition rate 6g→5g [1/yr].

    Returns
    -------
    float
        dN_e/dt [1/yr].
    """

    return gamma_sr_yr * Ne - gamma_t_yr * Ne * Ng


def dNgdt(Ne, Ng, gamma_sr_yr, gamma_t_yr):
    """
    Time derivative dN_e/dt in 1/yr units.

    Parameters
    ----------
    Ne : float
        Excited-state population N_e (dimensionless count).
    Ng : float
        Ground-state population N_g (dimensionless count).
    gamma_sr_yr : float
        Superradiance growth rate for the excited state [1/yr].
    gamma_t_yr : float
        Transition rate 6g→5g [1/yr].

    Returns
    -------
    float
        dN_e/dt [1/yr].
    """
    return gamma_sr_yr * Ng + gamma_t_yr * Ne * Ng


def strain_envelope(dist_GeV_inv, alpha, ne, ng, mu, gamma_t_GeV, Ng, Ne):
    """
    Strain envelope amplitude (no oscillatory factor).

    Parameters
    ----------
    dist_GeV_inv : float
        Source distance in natural units [GeV^-1] (e.g. kpc_to_GeV * d_kpc).
    alpha : float
        Dimensionless coupling α = G M μ.
    ne : int
        Excited-state principal quantum number.
    ng : int
        Ground-state principal quantum number.
    mu : float
        Boson mass μ [GeV].
    gamma_t_GeV : float
        Transition rate Γ_t [GeV].
    Ng : float
        Ground-state population N_g.
    Ne : float
        Excited-state population N_e.

    Returns
    -------
    float
        Strain envelope h(t) (dimensionless).
    """
    omega_tr = 0.5 * mu * alpha**2 * ((1.0 / ng**2) - (1.0 / ne**2))
    # print("\nomega_tr =", omega_tr)
    # print("gamma_t =", gamma_t_GeV)
    # print("Ng =", Ng, ", Ne =", Ne)
    # print("G =", G, ", dist =", dist_GeV_inv, "\n")
    amp = np.sqrt(4.0 * G / (dist_GeV_inv**2 * omega_tr) *
                  gamma_t_GeV * Ng * Ne)
    return amp


# ============================================================
# FFT and Lorentzian fit
# ============================================================
def fft_continuous_from_time(time_yr, h_t, n_fft=2**20):
    """
    Strain envelope amplitude (no oscillatory factor).

    Parameters
    ----------
    dist_GeV_inv : float
        Source distance in natural units [GeV^-1] (e.g. kpc_to_GeV * d_kpc).
    alpha : float
        Dimensionless coupling α = G M μ.
    ne : int
        Excited-state principal quantum number.
    ng : int
        Ground-state principal quantum number.
    mu : float
        Boson mass μ [GeV].
    gamma_t_GeV : float
        Transition rate Γ_t [GeV].
    Ng : float
        Ground-state population N_g.
    Ne : float
        Excited-state population N_e.

    Returns
    -------
    float
        Strain envelope h(t) (dimensionless).
    """
    time_yr = np.asarray(time_yr, dtype=float)
    h_t     = np.asarray(h_t,      dtype=complex)

    t_s = time_yr * S_PER_YR
    t_min, t_max = t_s[0], t_s[-1]
    t_uniform = np.linspace(t_min, t_max, n_fft)

    h_real = np.interp(t_uniform, t_s, h_t.real)
    h_imag = np.interp(t_uniform, t_s, h_t.imag)
    h_uniform = h_real + 1j*h_imag
    h_uniform -= np.mean(h_uniform)

    dt = t_uniform[1] - t_uniform[0]
    H_f  = np.fft.fft(h_uniform) * dt
    f_Hz = np.fft.fftfreq(n_fft, d=dt)
    return f_Hz, H_f

from scipy.optimize import curve_fit
import numpy as np

def lorentzian(f, A, f0, gamma, C):
    """
    L(f) = A / ((f - f0)^2 + gamma^2) + C
    """
    return A / ((f - f0)**2 + gamma**2) + C

def log_lorentzian(f, logA, f0, loggamma, logC):
    """
    Model for log10 |H(f)| so we can fit in log space.

    log10 L(f) = log10( A / ((f - f0)^2 + gamma^2) + C )

    where A = 10^logA, gamma = 10^loggamma, C = 10^logC.
    """
    A     = 10.0**logA
    gamma = 10.0**loggamma
    C     = 10.0**logC
    L = A / ((f - f0)**2 + gamma**2) + C
    L = np.maximum(L, 1e-300)
    return np.log10(L)

def fit_lorentzian_to_fft(time_yr, h_t, n_fft=2**20, n_top=400):
    """
    1. Compute FFT of h(t)
    2. Keep positive frequencies
    3. Select the n_top largest points in |H(f)| (around the main line)
    4. Fit a Lorentzian with floor in log space to those points
    5. Return full positive spectrum and a *floorless* Lorentzian
       evaluated on it (so the tails decay instead of flattening).

    Returns
    -------
    f_pos : ndarray
        Positive frequency grid [Hz].
    H_pos : ndarray
        FFT[h(t)] on f_pos (complex, continuous-normalized).
    L_full : ndarray
        Best-fit Lorentzian evaluated on f_pos, with floor removed:
            L_full(f) = A_fit / ((f - f0_fit)^2 + gamma_fit^2)
    lorentz_params : tuple
        (A_fit, f0_fit, gamma_fit, C_fit) from the underlying fit.
    """
    # 1) FFT
    f_Hz, H_f = fft_continuous_from_time(time_yr, h_t, n_fft=n_fft)

    # 2) positive frequencies only
    mask_pos  = f_Hz > 0
    f_pos     = f_Hz[mask_pos]
    H_pos     = H_f[mask_pos]
    abs_H_pos = np.abs(H_pos)

    # 3) pick n_top largest amplitudes
    idx_sorted = np.argsort(abs_H_pos)[::-1]
    n_top = min(n_top, len(idx_sorted))
    idx_fit = idx_sorted[:n_top]

    f_fit  = f_pos[idx_fit]
    H_fit  = abs_H_pos[idx_fit]

    # sort by frequency for stability
    sort_idx = np.argsort(f_fit)
    f_fit  = f_fit[sort_idx]
    H_fit  = H_fit[sort_idx]

    # 4) initial guesses
    f0_guess = f_fit[np.argmax(H_fit)]
    df       = np.median(np.diff(np.sort(f_pos)))
    A0       = (H_fit.max() - H_fit.min()) * (df**2)
    gamma0   = df * 5.0
    C0       = H_fit.min()

    logA0     = np.log10(max(A0,     1e-40))
    loggamma0 = np.log10(max(gamma0, 1e-40))
    logC0     = np.log10(max(C0,     1e-40))
    p0 = [logA0, f0_guess, loggamma0, logC0]

    # fit log10 |H(f)| vs f
    ydata = np.log10(np.maximum(H_fit, 1e-300))
    popt, pcov = curve_fit(log_lorentzian, f_fit, ydata, p0=p0, maxfev=20000)

    logA_fit, f0_fit, loggamma_fit, logC_fit = popt
    A_fit     = 10.0**logA_fit
    gamma_fit = 10.0**loggamma_fit
    C_fit     = 10.0**logC_fit

    # 5) evaluate *floorless* Lorentzian for plotting/usage
    L_full = A_fit / ((f_pos - f0_fit)**2 + gamma_fit**2)

    print(f"Fitted center frequency f0   = {f0_fit:.6e} Hz")
    print(f"Fitted width gamma           = {gamma_fit:.6e} Hz")
    print(f"Estimated FWHM ≃ 2*gamma     = {2*gamma_fit:.6e} Hz")
    print(f"Fitted constant floor C      = {C_fit:.3e}")

    lorentz_params = (A_fit, f0_fit, gamma_fit, C_fit)
    return f_pos, H_pos, L_full, lorentz_params


# ============================================================
# Annihilation helpers (frequency-domain analytic strain)
# ============================================================

def _mp_E1_vec(z):
    """
    Vectorized exponential integral E1(z) using mpmath.

    Parameters
    ----------
    z : array_like of complex
        Argument of the exponential integral.

    Returns
    -------
    ndarray of complex
        E1(z) evaluated elementwise.
    """
    z = np.asarray(z, dtype=complex)
    return np.vectorize(lambda zz: complex(mp.e1(zz)))(z)


def omega_ann(mua, alpha, n):
    """
    Annihilation line angular frequency (in GeV).

    Parameters
    ----------
    mua : float
        Boson mass μ_a [GeV].
    alpha : float
        Gravitational fine-structure α = G M μ_a (dimensionless).
    n : int
        Principal quantum number.

    Returns
    -------
    float
        ω_ann ≈ 2 μ_a (1 - α^2 / (2 n^2)) [GeV].
    """
    return 2.0 * mua * (1.0 - alpha**2 / (2.0 * n**2))


def gamma_ann(l, alpha, M_GeV):
    """
    Annihilation decay rate Γ_a in GeV, following the approximate
    Arvanitaki-style scaling you coded.

    Parameters
    ----------
    l : int
        Orbital angular momentum quantum number.
    alpha : float
        Gravitational fine-structure constant α = G M μ (dimensionless).
    M_GeV : float
        Black-hole mass [GeV].

    Returns
    -------
    float
        Γ_a [GeV].
    """
    if l == 1:
        p = 17
    else:
        p = 4 * l + 1

    r_g = G * M_GeV          # GeV^-1
    return G * 1e-10 / r_g**3 * (((alpha / l) * 0.5)**p +
                                 ((alpha / l) * 0.5)**(p + 1))


def h_ann(delta_f_Hz, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV):
    """
    Core annihilation line shape h̃(Δf) for one complex amplitude.

    Parameters
    ----------
    delta_f_Hz : array_like
        Frequency offset from line center [Hz].
    mua : float
        Boson mass μ_a [GeV].
    M_solar : float
        BH mass in solar masses.
    n, l : int
        Principal and orbital quantum numbers.
    alpha : float
        Gravitational fine-structure α = G M μ (dimensionless).
    r_kpc : float
        Distance to source [kpc].
    omega_a_GeV : float
        Annihilation line energy [GeV].

    Returns
    -------
    ndarray of complex
        h̃(Δf) in strain/Hz.
    """
    # Convert to natural units
    M_GeV = M_solar * M_sun
    r     = r_kpc * kpc_to_GeV   # GeV^-1

    # Decay rate [GeV]
    Gamma = gamma_ann(l, alpha, M_GeV)

    # Max cloud population (dimensionless, Arvanitaki approx)
    N_max = 10.0**76 * (M_solar / 10.0)**2

    # Energy offset ΔE = h Δf  [GeV]
    delta_E = np.asarray(delta_f_Hz, dtype=float) * Hz_to_GeV

    # Argument of the exponential integral (dimensionless)
    z = 1j * delta_E / (Gamma * N_max)

    # Prefactor (in GeV^-1; overall gives strain/GeV before multiplying by h)
    pref = 1.0 / (2.0 * np.pi) * np.sqrt(4.0 * G /
                                         (Gamma * r**2 * omega_a_GeV))

    Ei_part = _mp_E1_vec(z)
    Ei      = np.exp(1j * delta_E / (Gamma * N_max)) * Ei_part

    # h̃_E has units strain/GeV; multiply by h [GeV·s] -> strain/Hz
    h_tilde = pref * np.exp(z) * Ei * h_GeV_s

    return h_tilde


def h_pcr_ann(f_grid_Hz, mua, M_solar, n, l, alpha, r_kpc, iota, phase):
    """
    Frequency-domain annihilation strain h̃_plus, h̃_cross on a given
    detector frequency grid.

    Parameters
    ----------
    f_grid_Hz : array_like
        Detector frequency grid [Hz].
    mua : float
        Boson mass μ_a [GeV].
    M_solar : float
        BH mass in solar masses.
    n, l : int
        Principal and orbital quantum numbers.
    alpha : float
        Gravitational fine-structure α = G M μ_a (dimensionless).
    r_kpc : float
        Distance to source [kpc].
    iota : float
        Inclination angle [rad].
    phase : float
        Overall phase offset.

    Returns
    -------
    h_plus : ndarray of complex
        Plus polarization h̃_+(f) [strain/Hz].
    h_cross : ndarray of complex
        Cross polarization h̃_×(f) [strain/Hz].
    f_grid_Hz : ndarray
        The input frequency grid [Hz].
    f_line_Hz : float
        Annihilation line frequency f_ann [Hz].
    """
    # Line energy and frequency
    omega_a_GeV = omega_ann(mua, alpha, n)
    f_a_Hz      = omega_a_GeV * GeV_to_Hz

    f_grid_Hz = np.asarray(f_grid_Hz, dtype=float)

    # Frequency offsets (two symmetric contributions)
    delta_f_minus = f_grid_Hz - f_a_Hz       # around +f_a
    delta_f_plus  = -(f_grid_Hz - f_a_Hz)    # around -f_a

    # Two analytic pieces
    hmin = h_ann(delta_f_minus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)
    hp   = h_ann(delta_f_plus,  mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)

    # Polarizations
    c = np.cos(iota)
    h_plus  = (1.0 + c**2) / 4.0 * (np.exp( 1j * phase) * hmin +
                                    np.exp(-1j * phase) * hp)
    h_cross = (c / (2.0j))        * (np.exp( 1j * phase) * hmin -
                                     np.exp(-1j * phase) * hp)

    return h_plus, h_cross, f_grid_Hz, f_a_Hz



# ============================================================
# MASTER FUNCTION for Level Transition
# Below is a function of the gravitational wave envelope
# The oscillatory term serves to only shift the entire envelope to be 
# centered at the carrier frequency
# This is neglected since the frequency array can simply be shifted to whatever 
# the omega_tr variable is defined as (convert to Hz)
# ============================================================

def iso_gatom_level_tr_strain(
    M_solar      = 1e-6,       # BH mass in solar masses
    a_spin       = 0.999999,   # dimensionless BH spin
    alpha        = 1.0,        # gravitational fine-structure α = G M μ
    ne           = 6,
    ng           = 5,
    m            = None,       # default m = l = ng-1
    distance_kpc = 10.0,       # source distance [kpc]
    N_e0         = 1.0,
    N_g0         = 1.0,
    n_time       = 100000,     # number of time steps
    n_fft        = 2**20,      # FFT length
    n_top        = 400,        # number of FFT points used for Lorentzian fit
    verbose      = True,
    ):
    """
    Evolve the level populations and strain envelope, then compute
    the frequency-domain strain and fit a Lorentzian to its magnitude.

    Parameters
    ----------
    M_solar : float
        Black-hole mass in solar masses M_sun.
    a_spin : float
        Dimensionless spin parameter (0 <= a <= 1).
    alpha : float
        Gravitational fine-structure constant α = G M μ (dimensionless).
    ne : int
        Excited state principal quantum number (e.g. 6 for 6g).
    ng : int
        Ground state principal quantum number (e.g. 5 for 5g).
    m : int or None
        Azimuthal quantum number m; default m = ℓ = ng - 1.
    distance_kpc : float
        Source distance in kpc.
    N_e0 : float
        Initial excited-state population N_e(0).
    N_g0 : float
        Initial ground-state population N_g(0).
    n_time : int
        Number of time steps used in the evolution (log-spaced).
    n_fft : int
        FFT length for the frequency-domain transform.
    n_top : int
        Number of FFT points to use when fitting the Lorentzian.
    verbose : bool
        If True, print diagnostic information.

    Returns
    -------
    dict
        Dictionary with keys:
          'time_yr'        : time grid [yr]
          'h_t'            : strain envelope h(t) (dimensionless)
          'N_e'            : excited-state population history
          'N_g'            : ground-state population history
          'f_pos'          : positive frequency grid [Hz]
          'H_pos'          : FFT of h(t) on f_pos (complex)
          'L_full'         : best-fit Lorentzian (no floor) on f_pos
          'lorentz_params' : (A, f0, gamma, C) from the fit
          'omega_tr_GeV'   : transition angular frequency ω_tr [GeV]
          'Mu_a'           : boson mass μ_a [GeV]
    """

    # quantum numbers
    l = ng - 1
    if m is None:
        m = l

    # BH and coupling, exactly as in your standalone script
    M    = M_solar * M_sun
    mu_a = alpha / (G * M)

    # superradiant rates (GeV), same as script
    gamma_sre_GeV = super_gamma(ne, l, m, mu_a, M, a_spin) # Extra factor since this SR rate is off for some reason compared to the source notebook used? cant figure out where
    gamma_srg_GeV = super_gamma(ng, l, m, mu_a, M, a_spin)

    # convert to 1/yr
    gamma_sre_yr = gamma_sre_GeV * GeV_to_yrinv
    gamma_srg_yr = gamma_srg_GeV * GeV_to_yrinv

    # transition angular frequency (GeV), same as script
    omega_tr = 0.5 * mu_a * alpha**2 * ((1.0 / ng**2) - (1.0 / ne**2))

    # transition rate (GeV)
    print("\n alpha =", alpha)
    print(" omega_tr =", omega_tr)
    print(" M =", M)
    gamma_t_ne_GeV = gamma_t_6g_to_5g(alpha, omega_tr, M)
    gamma_t_ne_yr  = gamma_t_ne_GeV * GeV_to_yrinv

    # distance in GeV^-1
    dist_GeV_inv = distance_kpc * kpc_to_GeV

    # time stepping: change spacing if there are plotting issues
    dt_arr = np.logspace(-8, 1, num=n_time)  # years
    time_yr = []
    yr = 0.0

    # initialise populations and strain
    Ne = N_e0
    Ng = N_g0

    N_e_hist = []
    N_g_hist = []
    h_hist   = []

    omega_sr = mu_a

    # evolution loop, copied in structure from your standalone code
    for j in dt_arr:
        r = r_plus(M, a_spin)
        Omega_H = Omega_plus(a_spin, r)

        # same superradiance condition
        if m * Omega_H > omega_sr:
            # derivatives in 1/yr
            dNe = dNedt(Ne, Ng, gamma_sre_yr, gamma_t_ne_yr)
            dNg = dNgdt(Ne, Ng, gamma_srg_yr, gamma_t_ne_yr)

            Ne_next = Ne + dNe * j
            Ng_next = Ng + dNg * j

            if Ne_next < 0 or Ng_next < 0:
                break

            Ne, Ng = Ne_next, Ng_next

            yr += j
            #print("\n gamma_t_ne_GeV =", gamma_t_ne_GeV)

            h_val = strain_envelope(dist_GeV_inv, alpha, ne, ng,
                                    mu_a, gamma_t_ne_GeV, Ng, Ne)

            N_e_hist.append(Ne)
            N_g_hist.append(Ng)
            h_hist.append(h_val)
            time_yr.append(yr)
        else:
            # if SR condition fails, stop evolving (same logic)
            break

    time_yr = np.array(time_yr)
    N_e_hist = np.array(N_e_hist)
    N_g_hist = np.array(N_g_hist)
    h_hist   = np.array(h_hist, dtype=float)

# Verbose used for debugging
    if verbose:
        print(f"Max N_e ≈ {N_e_hist.max():.3e}")
        print(f"Max N_g ≈ {N_g_hist.max():.3e}")
        print(f"Max h(t) ≈ {np.max(np.abs(h_hist)):.3e}")
        print(f"omega_tr ≈ {omega_tr:.3e} GeV")
        print(f"gamma_sre_yr ≈ {gamma_sre_yr:.3e} 1/yr")
        print(f"gamma_srg_yr ≈ {gamma_srg_yr:.3e} 1/yr")
        print(f"gamma_t_yr  ≈ {gamma_t_ne_yr:.3e} 1/yr")
        print(f"gamma_sre_gev  ≈ {gamma_srg_GeV:.3e} 1/Gev")
        print(f"Mu_a  ≈ {mu_a:.3e} GeV")

    # FFT and Lorentzian fit, unchanged
    f_pos, H_pos, L_full, lorentz_params = fit_lorentzian_to_fft(
        time_yr, h_hist, n_fft=n_fft, n_top=n_top
    )

    return {
        "time_yr": time_yr,        # years
        "h_t": h_hist,             # strain envelope in time domain (dimensionless)
        "N_e": N_e_hist,           # excited-state population
        "N_g": N_g_hist,           # ground-state population
        "f_pos": f_pos,            # frequency grid [Hz]
        "H_pos": H_pos,            # FFT of h(t) on f_pos
        "L_full": L_full,          # best-fit Lorentzian (Use this for any sensitivity measurements since its nice and continuous)
        "lorentz_params": lorentz_params,  # (A, f0, gamma, C)
        "omega_tr_GeV": omega_tr,  # transition frequency [GeV]
        "Mu_a": mu_a,              # boson mass [GeV]
    }


# ============================================================
# MASTER FUNCTION for annihilation line strain (frequency-domain only)
# ============================================================

def iso_gatom_ann_strain(
    M_solar      = 3.1e-4,   # BH mass in solar masses
    mua          = 2e-16,    # boson mass [GeV]
    n            = 4,
    l            = None,
    alpha        = None,     # if None, compute α = G M μ
    distance_kpc = 1.0,
    iota         = 0.0,
    phase        = 0.0,
    f_min_Hz     = 1.0,
    f_max_Hz     = 1.0e12,
    n_f          = 50000,
    verbose      = True,
):
    """
    Analytic frequency-domain strain for annihilations, with no time-domain
    evolution. This mirrors the interface of the level-transition code,
    but directly returns h̃(f).

    Parameters
    ----------
    M_solar : float
        Black-hole mass in solar masses.
    mua : float
        Boson mass μ_a [GeV].
    n : int
        Principal quantum number (e.g. 4 for 4ℓ).
    l : int or None
        Orbital angular momentum quantum number; if None, l = n - 1.
    alpha : float or None
        Gravitational fine-structure constant α = G M μ. If None, it is
        computed from M_solar and mua.
    distance_kpc : float
        Source distance [kpc].
    iota : float
        Inclination angle [rad].
    phase : float
        Overall phase offset.
    f_min_Hz, f_max_Hz : float
        Frequency range [Hz] for the detector grid.
    n_f : int
        Number of frequency samples (log-spaced).
    verbose : bool
        If True, print diagnostics.

    Returns
    -------
    dict
        {
          "f_Hz"      : frequency grid [Hz],
          "h_plus"    : h̃_+(f) [strain/Hz],
          "h_cross"   : h̃_×(f) [strain/Hz],
          "h_c"       : characteristic strain h_c(f),
          "f_line_Hz" : annihilation line frequency f_ann [Hz],
          "alpha"     : gravitational fine-structure α,
          "mua_GeV"   : μ_a [GeV],
          "M_solar"   : BH mass in solar masses,
          "distance_kpc" : distance [kpc],
          "n"         : principal quantum number,
          "l"         : orbital quantum number,
          "iota"      : inclination angle [rad],
          "phase"     : phase offset [rad],
        }
    """
    if l is None:
        l = n - 1

    # If alpha not given, compute α = G M μ
    if alpha is None:
        M_GeV = M_solar * M_sun
        alpha = G * M_GeV * mua

    # Log-spaced detector frequency grid
    f_grid_Hz = np.logspace(np.log10(f_min_Hz), np.log10(f_max_Hz), n_f)

    # Compute frequency-domain annihilation strain
    h_plus, h_cross, f_grid_Hz, f_line_Hz = h_pcr_ann(
        f_grid_Hz, mua, M_solar, n, l, alpha, distance_kpc, iota, phase
    )

    # Characteristic strain
    h_c = 2.0 * f_grid_Hz * np.sqrt(np.abs(h_plus)**2 + np.abs(h_cross)**2)

    if verbose:
        print(f"alpha ≈ {alpha:.3e}")
        print(f"Annihilation line f_ann ≈ {f_line_Hz:.3e} Hz")
        print(f"max |h_plus| ≈ {np.max(np.abs(h_plus)):.3e}")
        print(f"max h_c ≈ {np.max(h_c):.3e}")

    return {
        "f_Hz": f_grid_Hz,
        "h_plus": h_plus,
        "h_cross": h_cross,
        "h_c": h_c,
        "f_line_Hz": f_line_Hz,
        "alpha": alpha,
        "mua_GeV": mua,
        "M_solar": M_solar,
        "distance_kpc": distance_kpc,
        "n": n,
        "l": l,
        "iota": iota,
        "phase": phase,
    }
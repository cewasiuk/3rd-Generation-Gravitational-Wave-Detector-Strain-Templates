"""
constants.py
------------
All physical constants, unit conversions, and every helper function
shared by more than one strain-template module.

Imported by: level_transition.py, annihilation.py
"""

import math as ma
import numpy as np
from scipy.optimize import curve_fit


# ============================================================
# PHYSICAL CONSTANTS & UNIT CONVERSIONS
# ============================================================
Mp          = 1.220890e19           # Planck mass [GeV]
G           = 1.0 / Mp**2          # Newton's constant [GeV^-2]
hbar_GeV_s  = 6.582119569e-25      # reduced Planck constant ħ [GeV·s]
h_GeV_s     = 4.1357e-24           # Planck constant h [GeV·s]  (E = h f)
s_per_yr    = 3.154e7              # seconds per year [s]
S_PER_YR    = s_per_yr             # alias

GeV_to_yrinv = (1.0 / hbar_GeV_s) * s_per_yr   # GeV  →  1/yr
Hz_to_GeV    = h_GeV_s                           # 1 Hz  →  GeV  (E = h f)
GeV_to_Hz    = 1.0 / h_GeV_s                     # GeV  →  Hz

M_sun       = 1.98847e30 * 5.60958885e26   # solar mass [GeV]
Mpc_to_GeV  = 1.56e38                      # 1 Mpc in natural units [GeV^-1]
kpc_to_GeV  = Mpc_to_GeV / 1000.0         # 1 kpc in natural units [GeV^-1]


# ============================================================
# DERIVED GEOMETRY HELPERS
# (used by both level_transition and annihilation)
# ============================================================

def r_g(M):
    """
    Gravitational radius r_g = G M  [GeV^-1].

    Parameters
    ----------
    M : float — BH mass [GeV]
    """
    return G * M


def r_plus(M, a):
    """
    Outer horizon radius r_+ = G M (1 + sqrt(1 - a^2))  [GeV^-1].

    Parameters
    ----------
    M : float — BH mass [GeV]
    a : float — dimensionless spin parameter
    """
    return r_g(M) * (1.0 + np.sqrt(1.0 - a**2))


def Omega_plus(a, r):
    """
    Horizon angular frequency Ω_+ = a / (2 r_+)  [GeV].

    Parameters
    ----------
    a : float — dimensionless spin
    r : float — horizon radius r_+ [GeV^-1]
    """
    return a / (2.0 * r)


# ============================================================
# SUPERRADIANCE HELPERS
# (used by both level_transition and annihilation)
# ============================================================

def omega_bound(mu, alpha, n):
    """
    Hydrogenic bound-state frequency  [GeV].

    ω_n ≈ μ (1 - α^2 / (2 n^2))
    """
    return mu * (1.0 - alpha**2 / (2.0 * n**2))


def glm(a, m, l, r, omega):
    """
    g_{lm} product factor entering the small-alpha superradiance rate.

    Parameters
    ----------
    a     : float — dimensionless spin
    m     : int   — azimuthal quantum number
    l     : int   — orbital quantum number
    r     : float — horizon radius r_+ [GeV^-1]
    omega : float — mode frequency [GeV]
    """
    factor = (a * m - 2.0 * r * omega)**2
    g = 1.0
    for k in range(1, l + 1):
        g *= k**2 * (1.0 - a**2) + factor
    return g


def super_gamma_lowalpha(n, l, m, mu, M, a):
    """
    Small-α superradiance rate Γ_low  [GeV].

    Parameters
    ----------
    n, l, m : int   — principal, orbital, azimuthal quantum numbers
    mu      : float — boson mass [GeV]
    M       : float — BH mass [GeV]
    a       : float — dimensionless spin
    """
    alpha   = G * M * mu
    r_p     = r_plus(M, a)
    Omega_p = Omega_plus(a, r_p)
    omega_n = omega_bound(mu, alpha, n)

    C = (
        2.0**(4*l + 1)
        * ma.factorial(n + l)
        / (n**(2*l + 4) * ma.factorial(n - l - 1))
        * (ma.factorial(l) / (ma.factorial(2*l) * ma.factorial(2*l + 1)))**2
    )
    g_lm = glm(a, m, l, r_p, omega_n)

    return (
        2.0 * r_p / M
        * C
        * g_lm
        * (m * Omega_p - omega_n)
        * alpha**(4*l + 5)
        / G / 4.0
    )


def gamma_wkb(mu, M, a, kappa=3.7):
    """
    WKB superradiance rate Γ_WKB  [GeV].

    Parameters
    ----------
    mu    : float — boson mass [GeV]
    M     : float — BH mass [GeV]
    a     : float — dimensionless spin
    kappa : float — WKB exponent prefactor (default 3.7)
    """
    alpha = G * M * mu
    return 1.0e-7 / r_g(M) * np.exp(-kappa * alpha)


def super_gamma(n, l, m, mu, M, a, alpha_low=0.1, alpha_high=15.0):
    """
    Matched superradiance rate: low-α analytic blended with WKB  [GeV].

    Parameters
    ----------
    n, l, m    : int   — quantum numbers
    mu         : float — boson mass [GeV]
    M          : float — BH mass [GeV]
    a          : float — dimensionless spin
    alpha_low  : float — blend lower boundary (default 0.1)
    alpha_high : float — blend upper boundary (default 15.0)
    """
    alpha = G * M * mu
    kappa = l + 0.5

    if alpha <= alpha_low:
        return super_gamma_lowalpha(n, l, m, mu, M, a)
    if alpha >= alpha_high:
        return gamma_wkb(mu, M, a, kappa=kappa)

    gamma_low = super_gamma_lowalpha(n, l, m, mu, M, a)
    gamma_hig = gamma_wkb(mu, M, a, kappa=kappa)
    w = (alpha - alpha_low) / (alpha_high - alpha_low)
    return np.exp((1.0 - w) * np.log(gamma_low) + w * np.log(gamma_hig))


# ============================================================
# TRANSITION RATE  6g → 5g
# (used by level_transition; lives here because it depends only
#  on G and r_g, both of which are defined above)
# ============================================================

def gamma_t_6g_to_5g(alpha, omega_tr, M):
    """
    Gravitational-wave transition rate Γ_t for 6g → 5g (ℓ = m = 4)  [GeV].

    Parameters
    ----------
    alpha    : float — dimensionless coupling α = G M μ
    omega_tr : float — transition angular frequency [GeV]
    M        : float — BH mass [GeV]
    """
    C  = (2**28 * 3**4 * 5**5) / (11**22 * np.pi)
    C *= (32 * np.pi / 15.0)          # angular integral ∫ sin^4θ dΩ
    P_t = C * G * alpha**12 / r_g(M)**4
    return P_t / omega_tr


# ============================================================
# POPULATION ODEs
# (used by level_transition)
# ============================================================

def dNedt(Ne, Ng, gamma_sr_yr, gamma_t_yr):
    """
    dN_e/dt  [1/yr].

    Parameters
    ----------
    Ne, Ng      : float — excited / ground state populations
    gamma_sr_yr : float — superradiance growth rate for excited state [1/yr]
    gamma_t_yr  : float — 6g→5g transition rate [1/yr]
    """
    return gamma_sr_yr * Ne - gamma_t_yr * Ne * Ng


def dNgdt(Ne, Ng, gamma_sr_yr, gamma_t_yr):
    """
    dN_g/dt  [1/yr].

    Parameters
    ----------
    Ne, Ng      : float — excited / ground state populations
    gamma_sr_yr : float — superradiance growth rate for ground state [1/yr]
    gamma_t_yr  : float — 6g→5g transition rate [1/yr]
    """
    return gamma_sr_yr * Ng + gamma_t_yr * Ne * Ng


# ============================================================
# STRAIN ENVELOPE
# (used by level_transition)
# ============================================================

def strain_envelope(dist_GeV_inv, alpha, ne, ng, mu, gamma_t_GeV, Ng, Ne):
    """
    Time-domain strain envelope amplitude h(t) (dimensionless).

    Parameters
    ----------
    dist_GeV_inv : float — source distance [GeV^-1]
    alpha        : float — dimensionless coupling α = G M μ
    ne, ng       : int   — excited / ground principal quantum numbers
    mu           : float — boson mass [GeV]
    gamma_t_GeV  : float — transition rate Γ_t [GeV]
    Ng, Ne       : float — ground / excited state populations
    """
    omega_tr = 0.5 * mu * alpha**2 * ((1.0 / ng**2) - (1.0 / ne**2))
    return np.sqrt(
        4.0 * G / (dist_GeV_inv**2 * omega_tr) * gamma_t_GeV * Ng * Ne
    )


# ============================================================
# FFT UTILITY
# (used by level_transition)
# ============================================================

def fft_continuous_from_time(time_yr, h_t, n_fft=2**20):
    """
    Interpolate h(t) onto a uniform grid and return its continuous FFT.

    Parameters
    ----------
    time_yr : array_like — time values [yr]
    h_t     : array_like — strain values (real or complex)
    n_fft   : int        — FFT length

    Returns
    -------
    f_Hz : ndarray — frequency grid [Hz]
    H_f  : ndarray — continuous-normalised FFT (complex)
    """
    time_yr   = np.asarray(time_yr, dtype=float)
    h_t       = np.asarray(h_t,     dtype=complex)
    t_s       = time_yr * S_PER_YR
    t_uniform = np.linspace(t_s[0], t_s[-1], n_fft)

    h_uniform  = np.interp(t_uniform, t_s, h_t.real).astype(complex)
    h_uniform += 1j * np.interp(t_uniform, t_s, h_t.imag)
    h_uniform -= np.mean(h_uniform)

    dt   = t_uniform[1] - t_uniform[0]
    H_f  = np.fft.fft(h_uniform) * dt
    f_Hz = np.fft.fftfreq(n_fft, d=dt)
    return f_Hz, H_f


# ============================================================
# LORENTZIAN FITTING
# (used by level_transition)
# ============================================================

def lorentzian(f, A, f0, gamma, C):
    """Standard Lorentzian with floor: A / ((f-f0)^2 + gamma^2) + C."""
    return A / ((f - f0)**2 + gamma**2) + C


def log_lorentzian(f, logA, f0, loggamma, logC):
    """log10 of a Lorentzian with floor, for stable log-space fitting."""
    A     = 10.0**logA
    gamma = 10.0**loggamma
    C     = 10.0**logC
    return np.log10(np.maximum(A / ((f - f0)**2 + gamma**2) + C, 1e-300))


def fit_lorentzian_to_fft(time_yr, h_t, n_fft=2**20, n_top=400):
    """
    FFT h(t), select the n_top largest spectral points, fit a Lorentzian
    in log space, and return the floorless best-fit curve.

    Parameters
    ----------
    time_yr : array_like — time grid [yr]
    h_t     : array_like — strain envelope (real or complex)
    n_fft   : int        — FFT length
    n_top   : int        — number of peak points used for the fit

    Returns
    -------
    f_pos          : ndarray — positive frequency grid [Hz]
    H_pos          : ndarray — FFT on f_pos (complex)
    L_full         : ndarray — best-fit floorless Lorentzian on f_pos
    lorentz_params : tuple   — (A, f0, gamma, C)
    """
    f_Hz, H_f = fft_continuous_from_time(time_yr, h_t, n_fft=n_fft)

    mask_pos  = f_Hz > 0
    f_pos     = f_Hz[mask_pos]
    H_pos     = H_f[mask_pos]
    abs_H_pos = np.abs(H_pos)

    idx_sorted = np.argsort(abs_H_pos)[::-1]
    n_top      = min(n_top, len(idx_sorted))
    f_fit      = f_pos[idx_sorted[:n_top]]
    H_fit      = abs_H_pos[idx_sorted[:n_top]]

    sort_idx     = np.argsort(f_fit)
    f_fit, H_fit = f_fit[sort_idx], H_fit[sort_idx]

    f0_guess = f_fit[np.argmax(H_fit)]
    df       = np.median(np.diff(np.sort(f_pos)))
    A0       = (H_fit.max() - H_fit.min()) * df**2
    gamma0   = df * 5.0
    C0       = H_fit.min()

    p0    = [np.log10(max(A0, 1e-40)), f0_guess,
             np.log10(max(gamma0, 1e-40)), np.log10(max(C0, 1e-40))]
    ydata = np.log10(np.maximum(H_fit, 1e-300))
    popt, _ = curve_fit(log_lorentzian, f_fit, ydata, p0=p0, maxfev=20000)

    logA_fit, f0_fit, loggamma_fit, logC_fit = popt
    A_fit     = 10.0**logA_fit
    gamma_fit = 10.0**loggamma_fit
    C_fit     = 10.0**logC_fit

    L_full = A_fit / ((f_pos - f0_fit)**2 + gamma_fit**2)

    print(f"Fitted centre frequency f0 = {f0_fit:.6e} Hz")
    print(f"Fitted width gamma         = {gamma_fit:.6e} Hz")
    print(f"Estimated FWHM ≃ 2·gamma   = {2*gamma_fit:.6e} Hz")
    print(f"Fitted floor C             = {C_fit:.3e}")

    return f_pos, H_pos, L_full, (A_fit, f0_fit, gamma_fit, C_fit)
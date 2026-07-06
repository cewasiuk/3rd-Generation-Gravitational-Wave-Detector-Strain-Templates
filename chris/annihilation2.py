

import warnings
import numpy as np

from constants import (
    G, M_sun, kpc_to_GeV,
    GeV_to_Hz, Hz_to_GeV, h_GeV_s,
    r_g,
)


HBAR_GeV_s = h_GeV_s / (2.0 * np.pi)          # reduced Planck constant [GeV*s]
AU_to_GeV  = kpc_to_GeV / 2.0626480624709636e8  # 1 kpc = 2.0626e8 AU



def alpha_of(M_solar, mua):
    """Gravitational fine-structure constant  alpha = G M mu  (dimensionless)."""
    return G * (M_solar * M_sun) * mua


def omega_GW(mua, alpha, n=2):

    return 2.0 * mua * (1.0 - alpha**2 / (2.0 * n**2))


def f_GW_Hz(mua, alpha, n=2):
    """Annihilation-line frequency [Hz]."""
    return omega_GW(mua, alpha, n) * GeV_to_Hz

# CLOUD MASS  (superradiance saturation)


CLOUD_EFFICIENCY_DEFAULT = 0.5   # -> M_a/M_b = 0.05 at alpha = 0.1

def cloud_mass(M_solar, alpha, cloud_efficiency=CLOUD_EFFICIENCY_DEFAULT):

    return cloud_efficiency * alpha * (M_solar * M_sun)


def superradiant_rate_211_leading(mua, alpha):

    return (1.0 / 24.0) * alpha**8 * mua



#  GW RADIATION POWER  (Yang & Huang Eq. 54)

def _P_num_J0(a):
    return (9.671e41 + 5.577e42 * a**2 + 1.474e43 * a**4 + 2.361e43 * a**6)

def _P_den(a):
    return (2.0 + a**2)**11 * (4.0 + a**2)**4


def P_ann(alpha, M_solar, M_a_GeV, a_star=0.0, _allow_spin=False):

    M_b_GeV = M_solar * M_sun
    a = alpha
    P_J0 = a**14 * (M_a_GeV / M_b_GeV)**2 * _P_num_J0(a) / _P_den(a)

    if a_star == 0.0:
        return P_J0

    if not _allow_spin:
        warnings.warn(
            "Test",
            RuntimeWarning,
        )
        return P_J0

    #: spin terms, J in geometric units J = a_star * G * M_b^2.
    # The relative magnitude of these vs the J=0 term has not been checked; do
    # not trust without reproducing Fig. 4 first.
    J = a_star * G * M_b_GeV**2
    B1 = a * (-3.839e80 - 2.111e81 * a**2 - 5.329e81 * a**4 - 8.165e81 * a**8)
    B2 = a**2 * (3.809e118 + 2.184e119 * a**2 + 5.799e119 * a**4 + 9.450e119 * a**6)
    pre = a**14 * (M_a_GeV / M_b_GeV)**2 / _P_den(a)
    # factor B0 back out of P_J0 to add sibling terms consistently
    term_J  = pre * (J / M_b_GeV**2) * B1
    term_J2 = pre * (J / M_b_GeV**2)**2 * B2
    return P_J0 + term_J + term_J2


def saturation_spin(alpha):

    x = 1.0 - alpha**2 / 8.0
    return 4.0 * alpha * x / (1.0 + 4.0 * alpha**2 * x**2)



# STRAIN  (power -> characteristic strain -> polarizations)

def h0_char(alpha, M_solar, M_a_GeV, distance_au, n=2, a_star=0.0):

    mua     = alpha / (G * M_solar * M_sun)        # invert alpha = G M mu
    omega_G = omega_GW(mua, alpha, n)              # GeV
    D       = distance_au * AU_to_GeV              # GeV^-1
    P       = P_ann(alpha, M_solar, M_a_GeV, a_star=a_star)
    return np.sqrt(10.0 * G * P) / (omega_G * D)


def polarization_amplitudes(h0, iota, phase=0.0):

    c = np.cos(iota)
    h_plus  = h0 * (1.0 + c**2) / 2.0 * np.exp(1j * phase)
    h_cross = h0 * c * 1j * np.exp(1j * phase)
    return h_plus, h_cross


#  TIMESCALES  (duration, drift, effective bandwidth)

def signal_duration(M_a_GeV, P_GeV2):

    return (M_a_GeV / P_GeV2) * HBAR_GeV_s


def freq_drift(alpha, M_solar, mua, n=2):

    # Calibration: anchor |fdot| at the benchmark and scale as alpha^17 * f_GW.
    f0 = f_GW_Hz(mua, alpha, n)

    return (alpha / 0.1)**17 * f0 * 1.0e-9   # Hz/s 


def effective_bandwidth(tau_s, fdot_Hz_s, T_obs_s):

    T_obs_eff = min(T_obs_s, tau_s)              # can't integrate longer than the signal
    return max(1.0 / T_obs_eff, abs(fdot_Hz_s) * T_obs_eff)



# SNR  (monochromatic, matched filter)

def snr_monochromatic(h_plus, h_cross, f_Hz, T_obs_s, S_n, tau_s=np.inf,
                      Fplus2=1.0, Fcross2=1.0):

    Sn = S_n(f_Hz) if callable(S_n) else S_n
    T_coh = min(T_obs_s, tau_s)
    power = (np.abs(h_plus)**2 * Fplus2 + np.abs(h_cross)**2 * Fcross2)
    return np.sqrt(power * T_coh / Sn)


#  RESOLVED-LINE SAMPLER  (for grid-based Fisher pipelines)
#
def sampled_strain(f_grid_Hz, f_GW, h_plus, h_cross, T_obs_s, tau_s=np.inf):

    f_grid_Hz = np.asarray(f_grid_Hz, dtype=float)
    T_coh = min(T_obs_s, tau_s)
    df = np.median(np.diff(f_grid_Hz))
    if df > 0.5 / T_coh:
        warnings.warn(
            f"test",
            RuntimeWarning,
        )
    # window normalised so int |w|^2 df = T_coh  
    w = T_coh * np.sinc((f_grid_Hz - f_GW) * T_coh)
    return h_plus * w, h_cross * w


#  MAIN FUNCTION

def iso_gatom_ann_strain(
    M_solar      = 1.22e-6,
    mua          = None,
    alpha        = None,
    n            = 2,
    distance_au  = 1.0,
    iota         = 0.0,
    phase        = 0.0,
    cloud_efficiency = CLOUD_EFFICIENCY_DEFAULT,
    M_a_GeV      = None,
    a_star       = 0.0,
    T_obs_s      = 4.0 * 86400.0,
    verbose      = True,
):

    M_b_GeV = M_solar * M_sun
    if alpha is None:
        if mua is None:
            raise ValueError("Provide either `mua` or `alpha`.")
        alpha = G * M_b_GeV * mua
    else:
        mua = alpha / (G * M_b_GeV)

    if M_a_GeV is None:
        M_a_GeV = cloud_mass(M_solar, alpha, cloud_efficiency)

    f_line = f_GW_Hz(mua, alpha, n)
    P      = P_ann(alpha, M_solar, M_a_GeV, a_star=a_star)
    h0     = h0_char(alpha, M_solar, M_a_GeV, distance_au, n=n, a_star=a_star)
    h_plus, h_cross = polarization_amplitudes(h0, iota, phase)
    tau    = signal_duration(M_a_GeV, P)
    fdot   = freq_drift(alpha, M_solar, mua, n)
    df_eff = effective_bandwidth(tau, fdot, T_obs_s)

    if verbose:
        print(f"alpha            = {alpha:.4e}")
        print(f"mu_a             = {mua:.4e} GeV")
        print(f"M_a / M_b        = {M_a_GeV / M_b_GeV:.4e}")
        print(f"f_line           = {f_line:.4e} Hz")
        print(f"P (GW power)     = {P:.4e} GeV^2")
        print(f"h0 (char strain) = {h0:.4e}")
        print(f"tau (duration)   = {tau:.4e} s = {tau/86400:.3f} days")
        print(f"fdot (drift)     = {fdot:.4e} Hz/s")
        print(f"Delta_f_eff      = {df_eff:.4e} Hz")

    return {
        "alpha":        alpha,
        "mua_GeV":      mua,
        "M_solar":      M_solar,
        "M_a_GeV":      M_a_GeV,
        "M_a_over_M_b": M_a_GeV / M_b_GeV,
        "distance_au":  distance_au,
        "n":            n,
        "iota":         iota,
        "phase":        phase,
        "a_star":       a_star,
        "f_line_Hz":    f_line,
        "P_GeV2":       P,
        "h0":           h0,
        "h_plus":       h_plus,
        "h_cross":      h_cross,
        "tau_s":        tau,
        "fdot_Hz_s":    fdot,
        "df_eff_Hz":    df_eff,
        "T_obs_s":      T_obs_s,
    }



# TESTS  (benchmarks from the three papers)

def _run_self_tests():
    print("=" * 64)
    print("SELF-TESTS vs published benchmarks")
    print("=" * 64)

    # CAPP benchmark point
    M_solar, alpha, dist_au = 1.22e-6, 0.1, 1.0
    M_a = cloud_mass(M_solar, alpha)

    # (1) line frequency vs CAPP Eq. 1 (~5.3 GHz for the matched axion mass)
    mua = alpha / (G * M_solar * M_sun)
    f = f_GW_Hz(mua, alpha)
    print(f"[1] f_line = {f:.3e} Hz   (CAPP Eq.1 ~5.3e9 Hz)   "
          f"{'PASS' if 4e9 < f < 7e9 else 'FAIL'}")

    # (2) characteristic strain vs CAPP Eq. 3 (~1.0e-22)
    h0 = h0_char(alpha, M_solar, M_a, dist_au)
    print(f"[2] h0     = {h0:.3e}    (CAPP Eq.3 ~1.0e-22)     "
          f"{'PASS' if 0.5e-22 < h0 < 2e-22 else 'FAIL'}")

    # (3) duration vs CAPP Eq. 4 (~4.7 days)
    P = P_ann(alpha, M_solar, M_a)
    tau_d = signal_duration(M_a, P) / 86400.0
    print(f"[3] tau    = {tau_d:.2f} days  (CAPP Eq.4 ~4.7 d)      "
          f"{'PASS' if 1.0 < tau_d < 15.0 else 'FAIL'}")

    # (4) power scaling P ~ alpha^14  (Yang-Huang Eq. 54)
    def Pexp(a):
        return P_ann(a, M_solar, cloud_mass(M_solar, a) / a)  # hold M_a fixed
    e = np.log(Pexp(0.2) / Pexp(0.05)) / np.log(0.2 / 0.05)
    print(f"[4] P ~ alpha^{e:.2f}    (Yang-Huang Eq.54 -> 14)    "
          f"{'PASS' if 13.5 < e < 14.5 else 'FAIL'}")

    # (5) with calibrated M_a(alpha)=0.5 alpha M_b: h0 ~ alpha^7, tau ~ alpha^-15
    def h0_a(a): return h0_char(a, M_solar, cloud_mass(M_solar, a), dist_au)
    def tau_a(a):
        Ma = cloud_mass(M_solar, a)
        return signal_duration(Ma, P_ann(a, M_solar, Ma))
    eh = np.log(h0_a(0.2)/h0_a(0.05)) / np.log(0.2/0.05)
    et = np.log(tau_a(0.2)/tau_a(0.05)) / np.log(0.2/0.05)
    print(f"[5] h0 ~ alpha^{eh:.2f} (CAPP Eq.3 ->7)   "
          f"tau ~ alpha^{et:.2f} (CAPP Eq.4 ->-15)   "
          f"{'PASS' if (6.5<eh<7.5 and -15.5<et<-14.5) else 'FAIL'}")

    # (6) monochromatic vs sampled-line SNR agree when grid resolves df_eff
    T_obs = 4.0 * 86400.0
    res = iso_gatom_ann_strain(M_solar=M_solar, alpha=alpha, distance_au=dist_au,
                               iota=0.4, T_obs_s=T_obs, verbose=False)
    Sn = 1e-48  # flat PSD [1/Hz], arbitrary for the cross-check
    rho_mono = snr_monochromatic(res["h_plus"], res["h_cross"],
                                 res["f_line_Hz"], T_obs, Sn, tau_s=res["tau_s"])
    Tcoh = min(T_obs, res["tau_s"])
    fG = res["f_line_Hz"]
    fgrid = np.linspace(fG - 50.0/Tcoh, fG + 50.0/Tcoh, 200001)
    hp, hc = sampled_strain(fgrid, fG, res["h_plus"], res["h_cross"], T_obs, res["tau_s"])
    rho_samp = np.sqrt(np.trapezoid((np.abs(hp)**2 + np.abs(hc)**2) / Sn, fgrid))
    rel = abs(rho_mono - rho_samp) / rho_mono
    print(f"[6] SNR mono={rho_mono:.3e} sampled={rho_samp:.3e} (rel.diff {rel:.1%})  "
          f"{'PASS' if rel < 0.05 else 'FAIL'}")
    print("=" * 64)


if __name__ == "__main__":
    print("Master-function demo at the CAPP benchmark:\n")
    iso_gatom_ann_strain(M_solar=1.22e-6, alpha=0.1, distance_au=1.0, iota=0.4)
    print()
    _run_self_tests()
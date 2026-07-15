if !isdefined(@__MODULE__, :ChrisWF)
    include("constants.jl")
end
using .ChrisWF
using Printf
using Statistics: median

const HBAR_GeV_s = h_GeV_s / (2.0 * pi)
const AU_to_GeV = kpc_to_GeV / 2.0626480624709636e8
const CLOUD_EFFICIENCY_DEFAULT = 0.5

"""Gravitational fine-structure constant `alpha = G M mu`."""
function alpha_of(M_solar, mua)
    return G * (M_solar * M_sun) * mua
end

function omega_GW(mua, alpha, n=2)
    return 2.0 * mua * (1.0 - alpha^2 / (2.0 * n^2))
end

"""Annihilation-line frequency in Hz."""
function f_GW_Hz(mua, alpha, n=2)
    return omega_GW(mua, alpha, n) * GeV_to_Hz
end

function cloud_mass(M_solar, alpha, cloud_efficiency=CLOUD_EFFICIENCY_DEFAULT)
    return cloud_efficiency * alpha * (M_solar * M_sun)
end

function superradiant_rate_211_leading(mua, alpha)
    return (1.0 / 24.0) * alpha^8 * mua
end

function _P_num_J0(a)
    return 9.671e41 + 5.577e42 * a^2 + 1.474e43 * a^4 + 2.361e43 * a^6
end

function _P_den(a)
    return (2.0 + a^2)^11 * (4.0 + a^2)^4
end

function P_ann(alpha, M_solar, M_a_GeV, a_star=0.0; _allow_spin=false)
    M_b_GeV = M_solar * M_sun
    a = alpha
    P_J0 = a^14 * (M_a_GeV / M_b_GeV)^2 * _P_num_J0(a) / _P_den(a)

    if a_star == 0.0
        return P_J0
    end

    if !_allow_spin
        @warn "Test"
        return P_J0
    end

    # Spin terms, with J in geometric units: J = a_star * G * M_b^2.
    J = a_star * G * M_b_GeV^2
    B1 = a * (-3.839e80 - 2.111e81 * a^2 - 5.329e81 * a^4 - 8.165e81 * a^8)
    B2 = a^2 * (3.809e118 + 2.184e119 * a^2 + 5.799e119 * a^4 + 9.450e119 * a^6)
    pre = a^14 * (M_a_GeV / M_b_GeV)^2 / _P_den(a)
    term_J = pre * (J / M_b_GeV^2) * B1
    term_J2 = pre * (J / M_b_GeV^2)^2 * B2
    return P_J0 + term_J + term_J2
end

function saturation_spin(alpha)
    x = 1.0 - alpha^2 / 8.0
    return 4.0 * alpha * x / (1.0 + 4.0 * alpha^2 * x^2)
end

function h0_char(alpha, M_solar, M_a_GeV, distance_au, n=2, a_star=0.0)
    mua = alpha / (G * M_solar * M_sun)
    omega_G = omega_GW(mua, alpha, n)
    D = distance_au * AU_to_GeV
    P = P_ann(alpha, M_solar, M_a_GeV, a_star)
    return sqrt(10.0 * G * P) / (omega_G * D)
end

function polarization_amplitudes(h0, iota, phase=0.0)
    c = cos(iota)
    h_plus = h0 * (1.0 + c^2) / 2.0 * exp(im * phase)
    h_cross = h0 * c * im * exp(im * phase)
    return h_plus, h_cross
end

function signal_duration(M_a_GeV, P_GeV2)
    return (M_a_GeV / P_GeV2) * HBAR_GeV_s
end

function freq_drift(alpha, M_solar, mua, n=2)
    f0 = f_GW_Hz(mua, alpha, n)
    return (alpha / 0.1)^17 * f0 * 1.0e-9
end

function effective_bandwidth(tau_s, fdot_Hz_s, T_obs_s)
    T_obs_eff = min(T_obs_s, tau_s)
    return max(1.0 / T_obs_eff, abs(fdot_Hz_s) * T_obs_eff)
end

function snr_monochromatic(
    h_plus,
    h_cross,
    f_Hz,
    T_obs_s,
    S_n,
    tau_s=Inf;
    Fplus2=1.0,
    Fcross2=1.0,
)
    Sn = applicable(S_n, f_Hz) ? S_n(f_Hz) : S_n
    T_coh = min(T_obs_s, tau_s)
    power = abs(h_plus)^2 * Fplus2 + abs(h_cross)^2 * Fcross2
    return sqrt(power * T_coh / Sn)
end

function sampled_strain(f_grid_Hz, f_GW, h_plus, h_cross, T_obs_s, tau_s=Inf)
    f_grid_Hz = collect(f_grid_Hz)
    T_coh = min(T_obs_s, tau_s)
    df = median(diff(f_grid_Hz))
    if df > 0.5 / T_coh
        @warn "test"
    end

    # Julia's sinc, like numpy.sinc, is sin(pi*x)/(pi*x).
    w = T_coh .* sinc.((f_grid_Hz .- f_GW) .* T_coh)
    return h_plus .* w, h_cross .* w
end

function iso_gatom_ann_strain(;
    M_solar=1.22e-6,
    mua=nothing,
    alpha=nothing,
    n=2,
    distance_au=1.0,
    iota=0.0,
    phase=0.0,
    cloud_efficiency=CLOUD_EFFICIENCY_DEFAULT,
    M_a_GeV=nothing,
    a_star=0.0,
    T_obs_s=4.0 * 86400.0,
    verbose=true,
)
    M_b_GeV = M_solar * M_sun
    if alpha === nothing
        if mua === nothing
            throw(ArgumentError("Provide either `mua` or `alpha`."))
        end
        alpha = G * M_b_GeV * mua
    else
        mua = alpha / (G * M_b_GeV)
    end

    if M_a_GeV === nothing
        M_a_GeV = cloud_mass(M_solar, alpha, cloud_efficiency)
    end

    f_line = f_GW_Hz(mua, alpha, n)
    P = P_ann(alpha, M_solar, M_a_GeV, a_star)
    h0 = h0_char(alpha, M_solar, M_a_GeV, distance_au, n, a_star)
    h_plus, h_cross = polarization_amplitudes(h0, iota, phase)
    tau = signal_duration(M_a_GeV, P)
    fdot = freq_drift(alpha, M_solar, mua, n)
    df_eff = effective_bandwidth(tau, fdot, T_obs_s)

    if verbose
        @printf("alpha            = %.4e\n", alpha)
        @printf("mu_a             = %.4e GeV\n", mua)
        @printf("M_a / M_b        = %.4e\n", M_a_GeV / M_b_GeV)
        @printf("f_line           = %.4e Hz\n", f_line)
        @printf("P (GW power)     = %.4e GeV^2\n", P)
        @printf("h0 (char strain) = %.4e\n", h0)
        @printf("tau (duration)   = %.4e s = %.3f days\n", tau, tau / 86400.0)
        @printf("fdot (drift)     = %.4e Hz/s\n", fdot)
        @printf("Delta_f_eff      = %.4e Hz\n", df_eff)
    end

    return Dict(
        "alpha" => alpha,
        "mua_GeV" => mua,
        "M_solar" => M_solar,
        "M_a_GeV" => M_a_GeV,
        "M_a_over_M_b" => M_a_GeV / M_b_GeV,
        "distance_au" => distance_au,
        "n" => n,
        "iota" => iota,
        "phase" => phase,
        "a_star" => a_star,
        "f_line_Hz" => f_line,
        "P_GeV2" => P,
        "h0" => h0,
        "h_plus" => h_plus,
        "h_cross" => h_cross,
        "tau_s" => tau,
        "fdot_Hz_s" => fdot,
        "df_eff_Hz" => df_eff,
        "T_obs_s" => T_obs_s,
    )
end

function _run_self_tests()
    println("="^64)
    println("SELF-TESTS vs published benchmarks")
    println("="^64)

    M_solar, alpha, dist_au = 1.22e-6, 0.1, 1.0
    M_a = cloud_mass(M_solar, alpha)

    mua = alpha / (G * M_solar * M_sun)
    f = f_GW_Hz(mua, alpha)
    @printf("[1] f_line = %.3e Hz   (CAPP Eq.1 ~5.3e9 Hz)   %s\n",
        f, 4e9 < f < 7e9 ? "PASS" : "FAIL")

    h0 = h0_char(alpha, M_solar, M_a, dist_au)
    @printf("[2] h0     = %.3e    (CAPP Eq.3 ~1.0e-22)     %s\n",
        h0, 0.5e-22 < h0 < 2e-22 ? "PASS" : "FAIL")

    P = P_ann(alpha, M_solar, M_a)
    tau_d = signal_duration(M_a, P) / 86400.0
    @printf("[3] tau    = %.2f days  (CAPP Eq.4 ~4.7 d)      %s\n",
        tau_d, 1.0 < tau_d < 15.0 ? "PASS" : "FAIL")

    Pexp(a) = P_ann(a, M_solar, cloud_mass(M_solar, a) / a)
    e = log(Pexp(0.2) / Pexp(0.05)) / log(0.2 / 0.05)
    @printf("[4] P ~ alpha^%.2f    (Yang-Huang Eq.54 -> 14)    %s\n",
        e, 13.5 < e < 14.5 ? "PASS" : "FAIL")

    h0_a(a) = h0_char(a, M_solar, cloud_mass(M_solar, a), dist_au)
    function tau_a(a)
        Ma = cloud_mass(M_solar, a)
        return signal_duration(Ma, P_ann(a, M_solar, Ma))
    end
    eh = log(h0_a(0.2) / h0_a(0.05)) / log(0.2 / 0.05)
    et = log(tau_a(0.2) / tau_a(0.05)) / log(0.2 / 0.05)
    pass5 = 6.5 < eh < 7.5 && -15.5 < et < -14.5
    @printf("[5] h0 ~ alpha^%.2f (CAPP Eq.3 ->7)   tau ~ alpha^%.2f (CAPP Eq.4 ->-15)   %s\n",
        eh, et, pass5 ? "PASS" : "FAIL")

    T_obs = 4.0 * 86400.0
    result = iso_gatom_ann_strain(;
        M_solar=M_solar,
        alpha=alpha,
        distance_au=dist_au,
        iota=0.4,
        T_obs_s=T_obs,
        verbose=false,
    )
    Sn = 1e-48
    rho_mono = snr_monochromatic(
        result["h_plus"], result["h_cross"], result["f_line_Hz"], T_obs, Sn,
        result["tau_s"],
    )
    Tcoh = min(T_obs, result["tau_s"])
    fG = result["f_line_Hz"]
    fgrid = collect(range(fG - 50.0 / Tcoh, fG + 50.0 / Tcoh, length=200001))
    hp, hc = sampled_strain(
        fgrid, fG, result["h_plus"], result["h_cross"], T_obs, result["tau_s"],
    )
    integrand = (abs2.(hp) .+ abs2.(hc)) ./ Sn
    rho_samp = sqrt(sum((integrand[1:end-1] .+ integrand[2:end]) .* diff(fgrid)) / 2.0)
    rel = abs(rho_mono - rho_samp) / rho_mono
    @printf("[6] SNR mono=%.3e sampled=%.3e (rel.diff %.1f%%)  %s\n",
        rho_mono, rho_samp, 100.0 * rel, rel < 0.05 ? "PASS" : "FAIL")
    println("="^64)

    return nothing
end

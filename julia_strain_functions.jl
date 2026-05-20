using LinearAlgebra
using FFTW
using SpecialFunctions
using Statistics

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
G   = 1.0 / Mp^2                   # GeV^-2
hbar_GeV_s = 6.582119569e-25       # GeV·s
s_per_yr   = 3.154e7
GeV_to_yrinv = (1.0 / hbar_GeV_s) * s_per_yr

M_sun     = 1.98847e30 * 5.60958885e26  # GeV
Mpc_to_GeV = 1.56e38
kpc_to_GeV = Mpc_to_GeV / 1000.0
S_PER_YR   = s_per_yr

h_GeV_s    = 4.1357e-24            # h in GeV·s (E = h f)
Hz_to_GeV  = h_GeV_s
GeV_to_Hz  = 1.0 / h_GeV_s


# ============================================================
# Isolated Level Transition helpers (frequency-domain analytic strain)
# ============================================================
function r_plus(M, a)
    """Horizon radius r_+ = G M (1 + sqrt(1 - a^2))."""
    return G * M * (1.0 + sqrt(1.0 - a^2))   # GeV^-1
end

function Omega_plus(a, r)
    """Horizon angular frequency Ω_+ = a / (2 r_+)."""
    return a / (2.0 * r)                     # GeV
end

function omega_bound(mu, alpha, n)
    """Hydrogenic bound-state frequency approximation."""
    return mu * (1.0 - alpha^2 / (2.0 * n^2))
end

function glm(a, m, l, r, omega)
    """g_{lm} product factor for the small-alpha rate."""
    factor = (a * m - 2.0 * r * omega)^2
    g = 1.0
    for k in 1:l
        g *= (k^2 * (1.0 - a^2) + factor)
    end
    return g
end


# ============================================================
# superradiance rates
# ============================================================
function super_gamma_lowalpha(n, l, m, mu, M, a)
    """Small-alpha superradiance rate Γ_low in GeV."""
    alpha = G * M * mu
    r_p   = r_plus(M, a)
    Omega_p = Omega_plus(a, r_p)
    omega_n = omega_bound(mu, alpha, n)

    C = (
        2.0^(4 * l + 1)
        * factorial(n + l)
        / (n^(2 * l + 4) * factorial(n - l - 1))
        * (factorial(l) / (factorial(2 * l) * factorial(2 * l + 1)))^2
    )
    g_lm = glm(a, m, l, r_p, omega_n)

    gamma = (
        2.0 * r_p / M
        * C
        * g_lm
        * (m * Omega_p - omega_n)
        * (alpha^(4 * l + 5))
        / G / 4.0
    )
    return gamma
end

function gamma_wkb(mu, M, a; kappa=3.7)
    """WKB superradiance rate Γ_WKB in GeV."""
    alpha = G * M * mu
    r_g   = G * M
    pref  = 1.0e-7 / r_g
    return pref * exp(-kappa * alpha)
end

function super_gamma(n, l, m, mu, M, a; alpha_low=0.1, alpha_high=15.0)
    """Matched superradiance rate (low-alpha + WKB) in GeV."""
    alpha = G * M * mu
    kappa = (l + 0.5)

    if alpha <= alpha_low
        return super_gamma_lowalpha(n, l, m, mu, M, a)
    end
    if alpha >= alpha_high
        return gamma_wkb(mu, M, a; kappa=kappa)
    end

    gamma_low = super_gamma_lowalpha(n, l, m, mu, M, a)
    gamma_hig = gamma_wkb(mu, M, a; kappa=kappa)
    log_low   = log(gamma_low)
    log_hig   = log(gamma_hig)
    w = (alpha - alpha_low) / (alpha_high - alpha_low)
    return exp((1.0 - w) * log_low + w * log_hig)
end


# ============================================================
# transition rate 6g -> 5g
# ============================================================
function gamma_t_6g_to_5g(alpha, omega_tr, M)
    """
    Transition rate Γ_t for 6g -> 5g (l = m = 4) in GeV.
    alpha: dimensionless coupling G M mu
    omega_tr: transition angular frequency in GeV
    """
    @show alpha
    @show omega_tr
    @show M
    r_g = G * M
    C = (2. ^28 * 3. ^4 * 5. ^5) / (11. ^22 * pi)
    C *= (32 * pi / 15.0)   # ∫ sin^4θ dΩ
    P_t = C * (G * alpha^12) / (r_g^4)
    return P_t / omega_tr
end


# ============================================================
# population ODEs and strain envelope
# ============================================================
function dNedt(Ne, Ng, gamma_sr_yr, gamma_t_yr)
    """
    Time derivative dN_e/dt in 1/yr units.
    """
    return gamma_sr_yr * Ne - gamma_t_yr * Ne * Ng
end

function dNgdt(Ne, Ng, gamma_sr_yr, gamma_t_yr)
    """
    Time derivative dN_g/dt in 1/yr units.
    """
    return gamma_sr_yr * Ng + gamma_t_yr * Ne * Ng
end

function strain_envelope(dist_GeV_inv, alpha, ne, ng, mu, gamma_t_GeV, Ng, Ne)
    """
    Strain envelope amplitude (no oscillatory factor).
    """
    @show dist_GeV_inv
    @show alpha
    @show ne
    @show ng
    @show mu
    @show gamma_t_GeV
    @show Ng
    @show Ne
    omega_tr = 0.5 * mu * alpha^2 * ((1.0 / ng^2) - (1.0 / ne^2))
    amp = sqrt(4.0 * G / (dist_GeV_inv^2 * omega_tr) *
               gamma_t_GeV * Ng * Ne)
    return amp
end


# ============================================================
# FFT and Lorentzian fit
# ============================================================
function fft_continuous_from_time(time_yr, h_t; n_fft=2^20)
    """
    Compute continuous-normalized FFT of a complex time series.
    """
    time_yr = collect(float.(time_yr))
    h_t     = ComplexF64.(h_t)

    t_s = time_yr .* S_PER_YR
    t_min, t_max = t_s[1], t_s[end]
    t_uniform = range(t_min, t_max; length=n_fft)

    h_real = interp1(t_s, real.(h_t), t_uniform)
    h_imag = interp1(t_s, imag.(h_t), t_uniform)
    h_uniform = h_real .+ 1im .* h_imag
    h_uniform .-= mean(h_uniform)

    dt = t_uniform[2] - t_uniform[1]
    H_f  = fft(h_uniform) .* dt
    f_Hz = fftfreq(n_fft, 1.0 / dt)
    return f_Hz, H_f
end

function lorentzian(f, A, f0, gamma, C)
    """L(f) = A / ((f - f0)^2 + gamma^2) + C"""
    return A ./ ((f .- f0).^2 .+ gamma^2) .+ C
end

function log_lorentzian(f, logA, f0, loggamma, logC)
    """
    Model for log10 |H(f)| so we can fit in log space.
    """
    A     = 10.0^logA
    gamma = 10.0^loggamma
    C     = 10.0^logC
    L = A ./ ((f .- f0).^2 .+ gamma^2) .+ C
    L = max.(L, 1e-300)
    return log10.(L)
end

function fit_lorentzian_to_fft(time_yr, h_t; n_fft=2^20, n_top=400)
    """
    1. Compute FFT of h(t)
    2. Keep positive frequencies
    3. Select the n_top largest points in |H(f)|
    4. Fit a Lorentzian with floor in log space to those points
    5. Return full positive spectrum and a floorless Lorentzian
    """
    f_Hz, H_f = fft_continuous_from_time(time_yr, h_t; n_fft=n_fft)

    mask_pos  = f_Hz .> 0
    f_pos     = f_Hz[mask_pos]
    H_pos     = H_f[mask_pos]
    abs_H_pos = abs.(H_pos)

    idx_sorted = sortperm(abs_H_pos; rev=true)
    n_top = min(n_top, length(idx_sorted))
    idx_fit = idx_sorted[1:n_top]

    f_fit  = f_pos[idx_fit]
    H_fit  = abs_H_pos[idx_fit]

    sort_idx = sortperm(f_fit)
    f_fit  = f_fit[sort_idx]
    H_fit  = H_fit[sort_idx]

    f0_guess = f_fit[argmax(H_fit)]
    df       = median(diff(sort(f_pos)))
    A0       = (maximum(H_fit) - minimum(H_fit)) * (df^2)
    gamma0   = df * 5.0
    C0       = minimum(H_fit)

    logA0     = log10(max(A0, 1e-40))
    loggamma0 = log10(max(gamma0, 1e-40))
    logC0     = log10(max(C0, 1e-40))
    p0 = [logA0, f0_guess, loggamma0, logC0]

    ydata = log10.(max.(H_fit, 1e-300))
    popt = curve_fit_log_lorentzian(f_fit, ydata, p0)

    logA_fit, f0_fit, loggamma_fit, logC_fit = popt
    A_fit     = 10.0^logA_fit
    gamma_fit = 10.0^loggamma_fit
    C_fit     = 10.0^logC_fit

    L_full = A_fit ./ ((f_pos .- f0_fit).^2 .+ gamma_fit^2)

    println("Fitted center frequency f0   = $(f0_fit) Hz")
    println("Fitted width gamma           = $(gamma_fit) Hz")
    println("Estimated FWHM ≃ 2*gamma     = $(2 * gamma_fit) Hz")
    println("Fitted constant floor C      = $(C_fit)")

    lorentz_params = (A_fit, f0_fit, gamma_fit, C_fit)
    return f_pos, H_pos, L_full, lorentz_params
end


# ============================================================
# Annihilation helpers (frequency-domain analytic strain)
# ============================================================
function _mp_E1_vec(z)
    """
    Vectorized exponential integral E1(z) using SpecialFunctions.expint.
    """
    z = ComplexF64.(z)
    return expint.(z)
end

function omega_ann(mua, alpha, n)
    """
    Annihilation line angular frequency (in GeV).
    """
    return 2.0 * mua * (1.0 - alpha^2 / (2.0 * n^2))
end

function gamma_ann(l, alpha, M_GeV)
    """
    Annihilation decay rate Γ_a in GeV.
    """
    p = l == 1 ? 17 : 4 * l + 1
    r_g = G * M_GeV          # GeV^-1
    return G * 1e-10 / r_g^3 * (((alpha / l) * 0.5)^p +
                                ((alpha / l) * 0.5)^(p + 1))
end

function h_ann(delta_f_Hz, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)
    """
    Core annihilation line shape h̃(Δf) for one complex amplitude.
    """
    M_GeV = M_solar * M_sun
    r     = r_kpc * kpc_to_GeV   # GeV^-1

    Gamma = gamma_ann(l, alpha, M_GeV)
    N_max = 10.0^76 * (M_solar / 10.0)^2

    delta_E = collect(float.(delta_f_Hz)) .* Hz_to_GeV
    z = 1im .* delta_E ./ (Gamma * N_max)

    pref = 1.0 / (2.0 * pi) * sqrt(4.0 * G /
                                   (Gamma * r^2 * omega_a_GeV))

    Ei_part = _mp_E1_vec(z)
    Ei      = exp.(1im .* delta_E ./ (Gamma * N_max)) .* Ei_part

    h_tilde = pref .* exp.(z) .* Ei .* h_GeV_s
    return h_tilde
end

function h_pcr_ann(f_grid_Hz, mua, M_solar, n, l, alpha, r_kpc, iota, phase)
    """
    Frequency-domain annihilation strain h̃_plus, h̃_cross on a grid.
    """
    omega_a_GeV = omega_ann(mua, alpha, n)
    f_a_Hz      = omega_a_GeV * GeV_to_Hz

    f_grid_Hz = collect(float.(f_grid_Hz))

    delta_f_minus = f_grid_Hz .- f_a_Hz
    delta_f_plus  = .-(f_grid_Hz .- f_a_Hz)

    hmin = h_ann(delta_f_minus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)
    hp   = h_ann(delta_f_plus,  mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)

    c = cos(iota)
    h_plus  = (1.0 + c^2) / 4.0 .* (exp( 1im * phase) .* hmin .+
                                    exp(-1im * phase) .* hp)
    h_cross = (c / (2.0im))       .* (exp( 1im * phase) .* hmin .-
                                     exp(-1im * phase) .* hp)

    return h_plus, h_cross, f_grid_Hz, f_a_Hz
end


# ============================================================
# MASTER FUNCTION for Level Transition
# ============================================================
function iso_gatom_level_tr_strain(; M_solar=1e-6, a_spin=0.999999, alpha=1.0,
                                   ne=6, ng=5, m=nothing, distance_kpc=10.0,
                                   N_e0=1.0, N_g0=1.0, n_time=100000,
                                   n_fft=2^20, n_top=400, verbose=true)
    """
    Evolve the level populations and strain envelope, then compute
    the frequency-domain strain and fit a Lorentzian to its magnitude.
    """
    l = ng - 1
    if m === nothing
        m = l
    end

    M    = M_solar * M_sun
    mu_a = alpha / (G * M)

    gamma_sre_GeV = super_gamma(ne, l, m, mu_a, M, a_spin)
    gamma_srg_GeV = super_gamma(ng, l, m, mu_a, M, a_spin)

    gamma_sre_yr = gamma_sre_GeV * GeV_to_yrinv
    gamma_srg_yr = gamma_srg_GeV * GeV_to_yrinv

    omega_tr = 0.5 * mu_a * alpha^2 * ((1.0 / ng^2) - (1.0 / ne^2))

    gamma_t_ne_GeV = gamma_t_6g_to_5g(alpha, omega_tr, M)
    @show gamma_t_ne_GeV
    gamma_t_ne_yr  = gamma_t_ne_GeV * GeV_to_yrinv

    dist_GeV_inv = distance_kpc * kpc_to_GeV

    dt_arr = 10.0 .^ range(-8, 1; length=n_time)
    time_yr = Float64[]
    yr = 0.0

    Ne = N_e0
    Ng = N_g0

    N_e_hist = Float64[]
    N_g_hist = Float64[]
    h_hist   = Float64[]

    omega_sr = mu_a

    for j in dt_arr
        r = r_plus(M, a_spin)
        Omega_H = Omega_plus(a_spin, r)

        if m * Omega_H > omega_sr
            dNe = dNedt(Ne, Ng, gamma_sre_yr, gamma_t_ne_yr)
            dNg = dNgdt(Ne, Ng, gamma_srg_yr, gamma_t_ne_yr)

            Ne_next = Ne + dNe * j
            Ng_next = Ng + dNg * j

            if Ne_next < 0 || Ng_next < 0
                break
            end

            Ne, Ng = Ne_next, Ng_next

            yr += j
            h_val = strain_envelope(dist_GeV_inv, alpha, ne, ng,
                                    mu_a, gamma_t_ne_GeV, Ng, Ne)

            push!(N_e_hist, Ne)
            push!(N_g_hist, Ng)
            push!(h_hist, h_val)
            push!(time_yr, yr)
        else
            break
        end
    end

    if verbose
        println("Max N_e ≈ $(maximum(N_e_hist))")
        println("Max N_g ≈ $(maximum(N_g_hist))")
        println("Max h(t) ≈ $(maximum(abs.(h_hist)))")
        println("omega_tr ≈ $(omega_tr) GeV")
        println("gamma_sre_yr ≈ $(gamma_sre_yr) 1/yr")
        println("gamma_srg_yr ≈ $(gamma_srg_yr) 1/yr")
        println("gamma_t_yr  ≈ $(gamma_t_ne_yr) 1/yr")
        println("gamma_sre_gev  ≈ $(gamma_srg_GeV) 1/Gev")
        println("Mu_a  ≈ $(mu_a) GeV")
    end

    f_pos, H_pos, L_full, lorentz_params = fit_lorentzian_to_fft(
        time_yr, h_hist; n_fft=n_fft, n_top=n_top
    )

    return Dict(
        "time_yr" => time_yr,
        "h_t" => h_hist,
        "N_e" => N_e_hist,
        "N_g" => N_g_hist,
        "f_pos" => f_pos,
        "H_pos" => H_pos,
        "L_full" => L_full,
        "lorentz_params" => lorentz_params,
        "omega_tr_GeV" => omega_tr,
        "Mu_a" => mu_a,
    )
end


# ============================================================
# MASTER FUNCTION for annihilation line strain (frequency-domain only)
# ============================================================
function iso_gatom_ann_strain(; M_solar=3.1e-4, mua=2e-16, n=4, l=nothing,
                              alpha=nothing, distance_kpc=1.0, iota=0.0,
                              phase=0.0, f_min_Hz=1.0, f_max_Hz=1.0e12,
                              n_f=50000, verbose=true)
    """
    Analytic frequency-domain strain for annihilations.
    """
    if l === nothing
        l = n - 1
    end

    if alpha === nothing
        M_GeV = M_solar * M_sun
        alpha = G * M_GeV * mua
    end

    f_grid_Hz = 10.0 .^ range(log10(f_min_Hz), log10(f_max_Hz); length=n_f)

    h_plus, h_cross, f_grid_Hz, f_line_Hz = h_pcr_ann(
        f_grid_Hz, mua, M_solar, n, l, alpha, distance_kpc, iota, phase
    )

    h_c = 2.0 .* f_grid_Hz .* sqrt.(abs.(h_plus).^2 .+ abs.(h_cross).^2)

    if verbose
        println("alpha ≈ $(alpha)")
        println("Annihilation line f_ann ≈ $(f_line_Hz) Hz")
        println("max |h_plus| ≈ $(maximum(abs.(h_plus)))")
        println("max h_c ≈ $(maximum(h_c))")
    end

    return Dict(
        "f_Hz" => f_grid_Hz,
        "h_plus" => h_plus,
        "h_cross" => h_cross,
        "h_c" => h_c,
        "f_line_Hz" => f_line_Hz,
        "alpha" => alpha,
        "mua_GeV" => mua,
        "M_solar" => M_solar,
        "distance_kpc" => distance_kpc,
        "n" => n,
        "l" => l,
        "iota" => iota,
        "phase" => phase,
    )
end


# ============================================================
# Minimal helpers to replace Python dependencies
# ============================================================
function interp1(x, y, xq)
    # Simple linear interpolation for monotonic x
    x = collect(x)
    y = collect(y)
    xq = collect(xq)
    yq = similar(xq, Float64)
    n = length(x)

    for i in eachindex(xq)
        xi = xq[i]
        if xi <= x[1]
            yq[i] = y[1]
        elseif xi >= x[end]
            yq[i] = y[end]
        else
            k = searchsortedfirst(x, xi)
            x0, x1 = x[k - 1], x[k]
            y0, y1 = y[k - 1], y[k]
            t = (xi - x0) / (x1 - x0)
            yq[i] = y0 + t * (y1 - y0)
        end
    end
    return yq
end

function fftfreq(n, fs)
    # fs is sampling rate (1/dt)
    val = fs / n
    n2 = fld(n, 2)
    if iseven(n)
        freqs = vcat(0:n2-1, -n2:-1) .* val
    else
        freqs = vcat(0:n2, -n2:-1) .* val
    end
    return freqs
end

function curve_fit_log_lorentzian(f_fit, ydata, p0)
    # Lightweight nonlinear fit using gradient-free search
    # This is a placeholder to mirror scipy.optimize.curve_fit.
    # For best results, replace with LsqFit.jl or Optim.jl.
    best = copy(p0)
    best_cost = _log_lorentzian_cost(f_fit, ydata, best)

    step = [0.1, 0.01 * maximum(f_fit), 0.1, 0.1]
    for _ in 1:2000
        trial = best .+ (rand(4) .- 0.5) .* step
        cost = _log_lorentzian_cost(f_fit, ydata, trial)
        if cost < best_cost
            best = trial
            best_cost = cost
        end
        step .*= 0.999
    end

    return best
end

function _log_lorentzian_cost(f_fit, ydata, p)
    logA, f0, loggamma, logC = p
    ypred = log_lorentzian(f_fit, logA, f0, loggamma, logC)
    return sum((ydata .- ypred).^2)
end

using SpecialFunctions
using FFTW
using LsqFit
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
const Mp = 1.220890e19                    # GeV
const G = 1.0 / Mp^2                      # GeV^-2
const hbar_GeV_s = 6.582119569e-25        # GeV·s
const s_per_yr = 3.154e7
const GeV_to_yrinv = (1.0 / hbar_GeV_s) * s_per_yr

const M_sun = 1.98847e30 * 5.60958885e26  # GeV
const Mpc_to_GeV = 1.56e38
const kpc_to_GeV = Mpc_to_GeV / 1000.0
const S_PER_YR = s_per_yr

const h_GeV_s = 4.1357e-24                # h in GeV·s (E = h f)
const Hz_to_GeV = h_GeV_s   
const GeV_to_Hz = 1.0 / h_GeV_s


# ============================================================
# Isolated Level Transition helpers (frequency-domain analytic strain)
# ============================================================
"""
    r_plus(M, a)

Horizon radius r_+ = G M (1 + sqrt(1 - a^2)).
"""
function r_plus(M::Real, a::Real)
    return G * M * (1.0 + sqrt(1.0 - a^2))   # GeV^-1
end


"""
    Omega_plus(a, r)

Horizon angular frequency Ω_+ = a / (2 r_+).
"""
function Omega_plus(a::Real, r::Real)
    return a / (2.0 * r)                     # GeV
end


"""
    omega_bound(mu, alpha, n)

Hydrogenic bound-state frequency approximation.
"""
function omega_bound(mu::Real, alpha::Real, n::Int)
    return mu * (1.0 - alpha^2 / (2.0 * n^2))
end


"""
    glm(a, m, l, r, omega)

g_{lm} product factor for the small-alpha rate.
"""
function glm(a::Real, m::Int, l::Int, r::Real, omega::Real)
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
"""
    super_gamma_lowalpha(n, l, m, mu, M, a)

Small-alpha superradiance rate Γ_low in GeV.
"""
function super_gamma_lowalpha(n::Int, l::Int, m::Int, mu::Real, M::Real, a::Real)
    alpha = G * M * mu
    r_p = r_plus(M, a)
    Omega_p = Omega_plus(a, r_p)
    omega_n = omega_bound(mu, alpha, n)

    C = (
        2.0^(4*l + 1)
        * factorial(n + l)
        / (n^(2*l + 4) * factorial(n - l - 1))
        * (factorial(l) / (factorial(2*l) * factorial(2*l + 1)))^2
    )
    g_lm = glm(a, m, l, r_p, omega_n)

    gamma = (
        2.0 * r_p / M
        * C
        * g_lm
        * (m * Omega_p - omega_n)
        * (alpha^(4*l + 5))
        / G / 4.0
    )
    return gamma
end


"""
    gamma_wkb(mu, M, a; kappa=3.7)

WKB superradiance rate Γ_WKB in GeV.
"""
function gamma_wkb(mu::Real, M::Real, a::Real; kappa::Real=3.7)
    alpha = G * M * mu
    r_g = G * M
    pref = 1.0e-7 / r_g
    return pref * exp(-kappa * alpha)
end


"""
    super_gamma(n, l, m, mu, M, a; alpha_low=0.1, alpha_high=15.0)

Matched superradiance rate (low-alpha + WKB) in GeV.
"""
function super_gamma(n::Int, l::Int, m::Int, mu::Real, M::Real, a::Real;
                     alpha_low::Real=0.1, alpha_high::Real=15.0)
    alpha = G * M * mu
    kappa = l + 0.5

    if alpha <= alpha_low
        return super_gamma_lowalpha(n, l, m, mu, M, a)
    end
    if alpha >= alpha_high
        return gamma_wkb(mu, M, a, kappa=kappa)
    end

    gamma_low = super_gamma_lowalpha(n, l, m, mu, M, a)
    gamma_hig = gamma_wkb(mu, M, a, kappa=kappa)
    log_low = log(gamma_low)
    log_hig = log(gamma_hig)
    w = (alpha - alpha_low) / (alpha_high - alpha_low)
    return exp((1.0 - w)*log_low + w*log_hig)
end


# ============================================================
# transition rate 6g -> 5g
# ============================================================
"""
    gamma_t_6g_to_5g(alpha, omega_tr, M)

Transition rate Γ_t for 6g → 5g (ℓ = m = 4) in GeV.

# Arguments
- `alpha`: dimensionless coupling G M mu
- `omega_tr`: transition angular frequency in GeV
- `M`: black-hole mass in GeV
"""
function gamma_t_6g_to_5g(alpha::Real, omega_tr::Real, M::Real)
    r_g = G * M
    C = (2. ^28 * 3. ^4 * 5. ^5) / (11. ^22 * π) * (32. *π/15.0) # ∫ sin^4θ dΩ
    P_t = C * (G * alpha^12) / (r_g^4)
    return P_t / omega_tr
end


# ============================================================
# population ODEs and strain envelope
# ============================================================
"""
    dNedt(Ne, Ng, gamma_sr_yr, gamma_t_yr)

Time derivative dN_e/dt in 1/yr units.

# Arguments
- `Ne`: Excited-state population N_e (dimensionless count)
- `Ng`: Ground-state population N_g (dimensionless count)
- `gamma_sr_yr`: Superradiance growth rate for the excited state [1/yr]
- `gamma_t_yr`: Transition rate 6g→5g [1/yr]

# Returns
- `Float64`: dN_e/dt [1/yr]
"""
function dNedt(Ne::Real, Ng::Real, gamma_sr_yr::Real, gamma_t_yr::Real)
    return gamma_sr_yr * Ne - gamma_t_yr * Ne * Ng
end


"""
    dNgdt(Ne, Ng, gamma_sr_yr, gamma_t_yr)

Time derivative dN_g/dt in 1/yr units.

# Arguments
- `Ne`: Excited-state population N_e (dimensionless count)
- `Ng`: Ground-state population N_g (dimensionless count)
- `gamma_sr_yr`: Superradiance growth rate for the excited state [1/yr]
- `gamma_t_yr`: Transition rate 6g→5g [1/yr]

# Returns
- `Float64`: dN_g/dt [1/yr]
"""
function dNgdt(Ne::Real, Ng::Real, gamma_sr_yr::Real, gamma_t_yr::Real)
    return gamma_sr_yr * Ng + gamma_t_yr * Ne * Ng
end


"""
    strain_envelope(dist_GeV_inv, alpha, ne, ng, mu, gamma_t_GeV, Ng, Ne)

Strain envelope amplitude (no oscillatory factor).

# Arguments
- `dist_GeV_inv`: Source distance in natural units [GeV^-1]
- `alpha`: Dimensionless coupling α = G M μ
- `ne`: Excited-state principal quantum number
- `ng`: Ground-state principal quantum number
- `mu`: Boson mass μ [GeV]
- `gamma_t_GeV`: Transition rate Γ_t [GeV]
- `Ng`: Ground-state population N_g
- `Ne`: Excited-state population N_e

# Returns
- `Float64`: Strain envelope h(t) (dimensionless)
"""
function strain_envelope(dist_GeV_inv::Real, alpha::Real, ne::Int, ng::Int,
                        mu::Real, gamma_t_GeV::Real, Ng::Real, Ne::Real)
    omega_tr = 0.5 * mu * alpha^2 * ((1.0 / ng^2) - (1.0 / ne^2))
    # println("omega_tr = $(omega_tr)")
    # println("gamma_t = $(gamma_t_GeV)")
    # println("Ng = $(Ng), Ne = $(Ne)")
    # println("G = $(G), dist = $(dist_GeV_inv)")
    amp = sqrt(4.0 * G / (dist_GeV_inv^2 * omega_tr) *
               gamma_t_GeV * Ng * Ne)
    return amp
end


# ============================================================
# FFT and Lorentzian fit
# ============================================================
"""
    fft_continuous_from_time(time_yr, h_t; n_fft=2^20)

Compute FFT of time-domain strain with continuous normalization.

# Arguments
- `time_yr`: Time array [yr]
- `h_t`: Strain time series (complex)
- `n_fft`: FFT length (default: 2^20)

# Returns
- `f_Hz`: Frequency array [Hz]
- `H_f`: FFT of h(t) [strain·s]
"""
function fft_continuous_from_time(time_yr::AbstractVector, h_t::AbstractVector;
                                 n_fft::Int=2^20)
    time_yr = Float64.(time_yr)
    h_t = ComplexF64.(h_t)

    t_s = time_yr .* S_PER_YR
    t_min, t_max = extrema(t_s)
    t_uniform = range(t_min, t_max, length=n_fft)

    # Interpolate real and imaginary parts
    h_real = linear_interp(t_s, real.(h_t), collect(t_uniform))
    h_imag = linear_interp(t_s, imag.(h_t), collect(t_uniform))
    h_uniform = h_real .+ 1im .* h_imag
    h_uniform .-= mean(h_uniform)

    dt = step(t_uniform)
    H_f = fft(h_uniform) .* dt
    f_Hz = fftfreq(n_fft, 1.0/dt)
    
    return f_Hz, H_f
end


"""
    linear_interp(x, y, x_new)

Simple linear interpolation helper.
"""
function linear_interp(x::AbstractVector, y::AbstractVector, x_new::AbstractVector)
    y_new = similar(x_new, eltype(y))
    for i in eachindex(x_new)
        xi = x_new[i]
        if xi <= x[1]
            y_new[i] = y[1]
        elseif xi >= x[end]
            y_new[i] = y[end]
        else
            idx = searchsortedfirst(x, xi)
            if idx > length(x)
                y_new[i] = y[end]
            elseif x[idx] == xi
                y_new[i] = y[idx]
            else
                x1, x2 = x[idx-1], x[idx]
                y1, y2 = y[idx-1], y[idx]
                y_new[i] = y1 + (y2 - y1) * (xi - x1) / (x2 - x1)
            end
        end
    end
    return y_new
end


"""
    lorentzian(f, p)

L(f) = A / ((f - f0)^2 + gamma^2) + C

where p = [A, f0, gamma, C]
"""
function lorentzian(f::Real, p::AbstractVector)
    A, f0, gamma, C = p
    return A / ((f - f0)^2 + gamma^2) + C
end


"""
    log_lorentzian(f, p)

Model for log10 |H(f)| in log space.

log10 L(f) = log10( A / ((f - f0)^2 + gamma^2) + C )

where p = [logA, f0, loggamma, logC]
"""
function log_lorentzian(f::Real, p::AbstractVector)
    logA, f0, loggamma, logC = p
    A = 10.0^logA
    gamma = 10.0^loggamma
    C = 10.0^logC
    L = A / ((f - f0)^2 + gamma^2) + C
    L = max(L, 1e-300)
    return log10(L)
end


"""
    fit_lorentzian_to_fft(time_yr, h_t; n_fft=2^20, n_top=400)

Compute FFT of h(t), fit a Lorentzian to the main peak, and return the spectrum.

# Returns
- `f_pos`: Positive frequency grid [Hz]
- `H_pos`: FFT[h(t)] on f_pos (complex)
- `L_full`: Best-fit Lorentzian on f_pos (no floor)
- `lorentz_params`: (A, f0, gamma, C) from fit
"""
function fit_lorentzian_to_fft(time_yr::AbstractVector, h_t::AbstractVector;
                               n_fft::Int=2^20, n_top::Int=400)
    # 1) FFT
    f_Hz, H_f = fft_continuous_from_time(time_yr, h_t, n_fft=n_fft)

    # 2) positive frequencies only
    mask_pos = f_Hz .> 0
    f_pos = f_Hz[mask_pos]
    H_pos = H_f[mask_pos]
    abs_H_pos = abs.(H_pos)

    # 3) pick n_top largest amplitudes
    idx_sorted = sortperm(abs_H_pos, rev=true)
    n_top = min(n_top, length(idx_sorted))
    idx_fit = idx_sorted[1:n_top]

    f_fit = f_pos[idx_fit]
    H_fit = abs_H_pos[idx_fit]

    # sort by frequency for stability
    sort_idx = sortperm(f_fit)
    f_fit = f_fit[sort_idx]
    H_fit = H_fit[sort_idx]

    # 4) initial guesses
    f0_guess = f_fit[argmax(H_fit)]
    df = median(diff(sort(f_pos)))
    A0 = (maximum(H_fit) - minimum(H_fit)) * (df^2)
    gamma0 = df * 5.0
    C0 = minimum(H_fit)

    logA0 = log10(max(A0, 1e-40))
    loggamma0 = log10(max(gamma0, 1e-40))
    logC0 = log10(max(C0, 1e-40))
    p0 = [logA0, f0_guess, loggamma0, logC0]

    # fit log10 |H(f)| vs f
    ydata = log10.(max.(H_fit, 1e-300))
    
    # Define model function for LsqFit
    model(f, p) = [log_lorentzian(fi, p) for fi in f]
    
    fit = curve_fit(model, f_fit, ydata, p0)
    popt = fit.param

    logA_fit, f0_fit, loggamma_fit, logC_fit = popt
    A_fit = 10.0^logA_fit
    gamma_fit = 10.0^loggamma_fit
    C_fit = 10.0^logC_fit

    # 5) evaluate *floorless* Lorentzian for plotting/usage
    L_full = @. A_fit / ((f_pos - f0_fit)^2 + gamma_fit^2)

    println("Fitted center frequency f0   = $(f0_fit) Hz")
    println("Fitted width gamma           = $(gamma_fit) Hz")
    println("Estimated FWHM ≃ 2*gamma     = $(2*gamma_fit) Hz")
    println("Fitted constant floor C      = $(C_fit)")

    lorentz_params = (A_fit, f0_fit, gamma_fit, C_fit)
    return f_pos, H_pos, L_full, lorentz_params
end


# ============================================================
# Annihilation helpers (frequency-domain analytic strain)
# ============================================================

"""
    _E1_vec(z)

Vectorized exponential integral E1(z) using SpecialFunctions.

# Arguments
- `z`: Array of complex numbers

# Returns
- Array of E1(z) evaluated elementwise
"""
function _E1_vec(z::AbstractVector{<:Complex})
    return [expint(zi) for zi in z]
end


"""
    omega_ann(mua, alpha, n)

Annihilation line angular frequency (in GeV).

# Arguments
- `mua`: Boson mass μ_a [GeV]
- `alpha`: Gravitational fine-structure α = G M μ_a (dimensionless)
- `n`: Principal quantum number

# Returns
- `Float64`: ω_ann ≈ 2 μ_a (1 - α^2 / (2 n^2)) [GeV]
"""
function omega_ann(mua::Real, alpha::Real, n::Int)
    return 2.0 * mua * (1.0 - alpha^2 / (2.0 * n^2))
end


"""
    gamma_ann(l, alpha, M_GeV)

Annihilation decay rate Γ_a in GeV.

# Arguments
- `l`: Orbital angular momentum quantum number
- `alpha`: Gravitational fine-structure constant α = G M μ (dimensionless)
- `M_GeV`: Black-hole mass [GeV]

# Returns
- `Float64`: Γ_a [GeV]
"""
function gamma_ann(l::Int, alpha::Real, M_GeV::Real)
    if l == 1
        p = 17
    else
        p = 4 * l + 1
    end

    r_g = G * M_GeV
    return G * 1e-10 / r_g^3 * (((alpha / l) * 0.5)^p +
                                ((alpha / l) * 0.5)^(p + 1))
end


"""
    h_ann(delta_f_Hz, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)

Core annihilation line shape h̃(Δf) for one complex amplitude.

# Arguments
- `delta_f_Hz`: Frequency offset from line center [Hz]
- `mua`: Boson mass μ_a [GeV]
- `M_solar`: BH mass in solar masses
- `n`, `l`: Principal and orbital quantum numbers
- `alpha`: Gravitational fine-structure α = G M μ (dimensionless)
- `r_kpc`: Distance to source [kpc]
- `omega_a_GeV`: Annihilation line energy [GeV]

# Returns
- Array of complex: h̃(Δf) in strain/Hz
"""
function h_ann(delta_f_Hz::AbstractVector, mua::Real, M_solar::Real, 
               n::Int, l::Int, alpha::Real, r_kpc::Real, omega_a_GeV::Real)
    # Convert to natural units
    M_GeV = M_solar * M_sun
    r = r_kpc * kpc_to_GeV   # GeV^-1

    # Decay rate [GeV]
    Gamma = gamma_ann(l, alpha, M_GeV)

    # Max cloud population (dimensionless, Arvanitaki approx)
    N_max = 10.0^76 * (M_solar / 10.0)^2

    # Energy offset ΔE = h Δf  [GeV]
    delta_E = Float64.(delta_f_Hz) .* Hz_to_GeV

    # Argument of the exponential integral (dimensionless)
    z = @. 1im * delta_E / (Gamma * N_max)

    # Prefactor (in GeV^-1; overall gives strain/GeV before multiplying by h)
    pref = 1.0 / (2.0 * π) * sqrt(4.0 * G /
                                  (Gamma * r^2 * omega_a_GeV))

    Ei_part = _E1_vec(z)
    Ei = @. exp(1im * delta_E / (Gamma * N_max)) * Ei_part

    # h̃_E has units strain/GeV; multiply by h [GeV·s] -> strain/Hz
    h_tilde = @. pref * exp(z) * Ei * h_GeV_s

    return h_tilde
end


"""
    h_pcr_ann(f_grid_Hz, mua, M_solar, n, l, alpha, r_kpc, iota, phase)

Frequency-domain annihilation strain h̃_plus, h̃_cross on a given detector frequency grid.

# Arguments
- `f_grid_Hz`: Detector frequency grid [Hz]
- `mua`: Boson mass μ_a [GeV]
- `M_solar`: BH mass in solar masses
- `n`, `l`: Principal and orbital quantum numbers
- `alpha`: Gravitational fine-structure α = G M μ_a (dimensionless)
- `r_kpc`: Distance to source [kpc]
- `iota`: Inclination angle [rad]
- `phase`: Overall phase offset

# Returns
- `h_plus`: Plus polarization h̃_+(f) [strain/Hz]
- `h_cross`: Cross polarization h̃_×(f) [strain/Hz]
- `f_grid_Hz`: The input frequency grid [Hz]
- `f_line_Hz`: Annihilation line frequency f_ann [Hz]
"""
function h_pcr_ann(f_grid_Hz::AbstractVector, mua::Real, M_solar::Real,
                   n::Int, l::Int, alpha::Real, r_kpc::Real,
                   iota::Real, phase::Real)
    # Line energy and frequency
    omega_a_GeV = omega_ann(mua, alpha, n)
    f_a_Hz = omega_a_GeV * GeV_to_Hz

    f_grid_Hz = Float64.(f_grid_Hz)

    # Frequency offsets (two symmetric contributions)
    delta_f_minus = f_grid_Hz .- f_a_Hz       # around +f_a
    delta_f_plus = -(f_grid_Hz .- f_a_Hz)     # around -f_a

    # Two analytic pieces
    hmin = h_ann(delta_f_minus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)
    hp = h_ann(delta_f_plus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)

    # Polarizations
    c = cos(iota)
    h_plus = @. (1.0 + c^2) / 4.0 * (exp(1im * phase) * hmin +
                                     exp(-1im * phase) * hp)
    h_cross = @. (c / (2.0im)) * (exp(1im * phase) * hmin -
                                  exp(-1im * phase) * hp)

    return h_plus, h_cross, f_grid_Hz, f_a_Hz
end


# ============================================================
# MASTER FUNCTION for Level Transition
# ============================================================

"""
    iso_gatom_level_tr_strain(; kwargs...)

Evolve the level populations and strain envelope, then compute
the frequency-domain strain and fit a Lorentzian to its magnitude.

# Keyword Arguments
- `M_solar=1e-6`: Black-hole mass in solar masses
- `a_spin=0.999999`: Dimensionless BH spin
- `alpha=1.0`: Gravitational fine-structure α = G M μ
- `ne=6`: Excited state principal quantum number
- `ng=5`: Ground state principal quantum number
- `m=nothing`: Azimuthal quantum number (default: m = l = ng-1)
- `distance_kpc=10.0`: Source distance [kpc]
- `N_e0=1.0`: Initial excited-state population
- `N_g0=1.0`: Initial ground-state population
- `n_time=100000`: Number of time steps
- `n_fft=2^20`: FFT length
- `n_top=400`: Number of FFT points for Lorentzian fit
- `verbose=true`: Print diagnostic information

# Returns
Dictionary with keys:
- `time_yr`: time grid [yr]
- `h_t`: strain envelope h(t) (dimensionless)
- `N_e`: excited-state population history
- `N_g`: ground-state population history
- `f_pos`: positive frequency grid [Hz]
- `H_pos`: FFT of h(t) on f_pos (complex)
- `L_full`: best-fit Lorentzian (no floor) on f_pos
- `lorentz_params`: (A, f0, gamma, C) from the fit
- `omega_tr_GeV`: transition angular frequency ω_tr [GeV]
- `Mu_a`: boson mass μ_a [GeV]
"""
function iso_gatom_level_tr_strain(;
    M_solar::Real = 1e-6,
    a_spin::Real = 0.999999,
    alpha::Real = 1.0,
    ne::Int = 6,
    ng::Int = 5,
    m::Union{Int,Nothing} = nothing,
    distance_kpc::Real = 10.0,
    N_e0::Real = 1.0,
    N_g0::Real = 1.0,
    n_time::Int = 100000,
    n_fft::Int = 2^20,
    n_top::Int = 400,
    verbose::Bool = true
    )

    # quantum numbers
    l = ng - 1
    if isnothing(m)
        m = l
    end

    # BH and coupling
    M = M_solar * M_sun
    mu_a = alpha / (G * M)

    # superradiant rates (GeV)
    gamma_sre_GeV = super_gamma(ne, l, m, mu_a, M, a_spin)
    gamma_srg_GeV = super_gamma(ng, l, m, mu_a, M, a_spin)

    # convert to 1/yr
    gamma_sre_yr = gamma_sre_GeV * GeV_to_yrinv
    gamma_srg_yr = gamma_srg_GeV * GeV_to_yrinv

    # transition angular frequency (GeV)
    omega_tr = 0.5 * mu_a * alpha^2 * ((1.0 / ng^2) - (1.0 / ne^2))

    # transition rate (GeV)
    gamma_t_ne_GeV = gamma_t_6g_to_5g(alpha, omega_tr, M)
    gamma_t_ne_yr = gamma_t_ne_GeV * GeV_to_yrinv

    # distance in GeV^-1
    dist_GeV_inv = distance_kpc * kpc_to_GeV

    # time stepping
    dt_arr = 10 .^ range(-8, 1, length=n_time)  # years
    time_yr = Float64[]
    yr = 0.0

    # initialise populations and strain
    Ne = N_e0
    Ng = N_g0

    N_e_hist = Float64[]
    N_g_hist = Float64[]
    h_hist = Float64[]

    omega_sr = mu_a

    # evolution loop
    for j in dt_arr
        r = r_plus(M, a_spin)
        Omega_H = Omega_plus(a_spin, r)

        # superradiance condition
        if m * Omega_H > omega_sr
            # derivatives in 1/yr
            dNe = dNedt(Ne, Ng, gamma_sre_yr, gamma_t_ne_yr)
            dNg_val = dNgdt(Ne, Ng, gamma_srg_yr, gamma_t_ne_yr)

            Ne_next = Ne + dNe * j
            Ng_next = Ng + dNg_val * j

            if Ne_next < 0 || Ng_next < 0
                break
            end

            Ne, Ng = Ne_next, Ng_next

            yr += j
            #println("gamma_t_ne_GeV =", gamma_t_ne_GeV)
            h_val = strain_envelope(dist_GeV_inv, alpha, ne, ng,
                                   mu_a, gamma_t_ne_GeV, Ng, Ne)

            push!(N_e_hist, Ne)
            push!(N_g_hist, Ng)
            push!(h_hist, h_val)
            push!(time_yr, yr)
        else
            # if SR condition fails, stop evolving
            break
        end
    end

    # Verbose for debugging
    if verbose
        println("Max N_e ≈ $(maximum(N_e_hist))")
        println("Max N_g ≈ $(maximum(N_g_hist))")
        println("Max h(t) ≈ $(maximum(abs.(h_hist)))")
        println("omega_tr ≈ $(omega_tr) GeV")
        println("gamma_sre_yr ≈ $(gamma_sre_yr) 1/yr")
        println("gamma_srg_yr ≈ $(gamma_srg_yr) 1/yr")
        println("gamma_t_yr  ≈ $(gamma_t_ne_yr) 1/yr")
        println("gamma_sre_gev  ≈ $(gamma_srg_GeV) 1/GeV")
        println("Mu_a  ≈ $(mu_a) GeV")
    end

    # FFT and Lorentzian fit
    f_pos, H_pos, L_full, lorentz_params = fit_lorentzian_to_fft(
        time_yr, h_hist, n_fft=n_fft, n_top=n_top
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
        "Mu_a" => mu_a
    )
end


# ============================================================
# MASTER FUNCTION for annihilation line strain (frequency-domain only)
# ============================================================

"""
    iso_gatom_ann_strain(; kwargs...)

Analytic frequency-domain strain for annihilations, with no time-domain evolution.

# Keyword Arguments
- `M_solar=3.1e-4`: Black-hole mass in solar masses
- `mua=2e-16`: Boson mass μ_a [GeV]
- `n=4`: Principal quantum number
- `l=nothing`: Orbital angular momentum (default: n-1)
- `alpha=nothing`: Gravitational fine-structure α (default: G M μ)
- `distance_kpc=1.0`: Source distance [kpc]
- `iota=0.0`: Inclination angle [rad]
- `phase=0.0`: Overall phase offset
- `f_min_Hz=1.0`: Minimum frequency [Hz]
- `f_max_Hz=1.0e12`: Maximum frequency [Hz]
- `n_f=50000`: Number of frequency samples
- `verbose=true`: Print diagnostics

# Returns
Dictionary with keys:
- `f_Hz`: frequency grid [Hz]
- `h_plus`: h̃_+(f) [strain/Hz]
- `h_cross`: h̃_×(f) [strain/Hz]
- `h_c`: characteristic strain h_c(f)
- `f_line_Hz`: annihilation line frequency f_ann [Hz]
- `alpha`: gravitational fine-structure α
- `mua_GeV`: μ_a [GeV]
- `M_solar`: BH mass in solar masses
- `distance_kpc`: distance [kpc]
- `n`: principal quantum number
- `l`: orbital quantum number
- `iota`: inclination angle [rad]
- `phase`: phase offset [rad]
"""
function iso_gatom_ann_strain(;
    M_solar::Real = 3.1e-4,
    mua::Real = 2e-16,
    n::Int = 4,
    l::Union{Int,Nothing} = nothing,
    alpha::Union{Real,Nothing} = nothing,
    distance_kpc::Real = 1.0,
    iota::Real = 0.0,
    phase::Real = 0.0,
    f_min_Hz::Real = 1.0,
    f_max_Hz::Real = 1.0e12,
    n_f::Int = 50000,
    verbose::Bool = true
    )

    if isnothing(l)
        l = n - 1
    end

    # If alpha not given, compute α = G M μ
    if isnothing(alpha)
        M_GeV = M_solar * M_sun
        alpha = G * M_GeV * mua
    end

    # Log-spaced detector frequency grid
    f_grid_Hz = 10 .^ range(log10(f_min_Hz), log10(f_max_Hz), length=n_f)

    # Compute frequency-domain annihilation strain
    h_plus, h_cross, f_grid_Hz, f_line_Hz = h_pcr_ann(
        f_grid_Hz, mua, M_solar, n, l, alpha, distance_kpc, iota, phase
    )

    # Characteristic strain
    h_c = @. 2.0 * f_grid_Hz * sqrt(abs(h_plus)^2 + abs(h_cross)^2)

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
        "phase" => phase
    )
end

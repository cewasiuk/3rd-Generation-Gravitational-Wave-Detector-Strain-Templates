module ChrisWF

using FFTW
using LinearAlgebra
using Printf
using SpecialFunctions
using Statistics

export Mp, G, hbar_GeV_s, h_GeV_s, s_per_yr, S_PER_YR, GeV_to_yrinv, Hz_to_GeV,
       GeV_to_Hz, M_sun, Mpc_to_GeV, kpc_to_GeV, r_g, r_plus, Omega_plus,
       omega_bound, glm, super_gamma_lowalpha, gamma_wkb, super_gamma,
       gamma_t_6g_to_5g, dNedt, dNgdt, strain_envelope,
    fft_continuous_from_time, lorentzian, log_lorentzian,
    fit_lorentzian_to_fft

const Mp           = 1.220890e19
const G            = 1.0 / Mp^2
const hbar_GeV_s   = 6.582119569e-25
const h_GeV_s      = 4.1357e-24
const s_per_yr     = 3.154e7
const S_PER_YR     = s_per_yr
const GeV_to_yrinv  = (1.0 / hbar_GeV_s) * s_per_yr
const Hz_to_GeV     = h_GeV_s
const GeV_to_Hz     = 1.0 / h_GeV_s
const M_sun        = 1.98847e30 * 5.60958885e26
const Mpc_to_GeV   = 1.56e38
const kpc_to_GeV   = Mpc_to_GeV / 1000.0

function r_g(M)
    return G * M
end

function r_plus(M, a)
    return r_g(M) * (1.0 + sqrt(1.0 - a^2))
end

function Omega_plus(a, r)
    return a / (2.0 * r)
end

function omega_bound(mu, alpha, n)
    return mu * (1.0 - alpha^2 / (2.0 * n^2))
end

function glm(a, m, l, r, omega)
    factor = (a * m - 2.0 * r * omega)^2
    g = 1.0
    for k in 1:l
        g *= k^2 * (1.0 - a^2) + factor
    end
    return g
end

function super_gamma_lowalpha(n, l, m, mu, M, a)
    alpha = G * M * mu
    r_p = r_plus(M, a)
    Omega_p = Omega_plus(a, r_p)
    omega_n = omega_bound(mu, alpha, n)

    C = (
        2.0^(4 * l + 1) * factorial(n + l) /
        (n^(2 * l + 4) * factorial(n - l - 1)) *
        (factorial(l) / (factorial(2 * l) * factorial(2 * l + 1)))^2
    )
    g_lm = glm(a, m, l, r_p, omega_n)

    return (
        2.0 * r_p / M * C * g_lm * (m * Omega_p - omega_n) * alpha^(4 * l + 5) / G / 4.0
    )
end

function gamma_wkb(mu, M, a; kappa=3.7)
    alpha = G * M * mu
    return 1.0e-7 / r_g(M) * exp(-kappa * alpha)
end

function super_gamma(n, l, m, mu, M, a; alpha_low=0.1, alpha_high=15.0)
    alpha = G * M * mu
    kappa = l + 0.5

    if alpha <= alpha_low
        return super_gamma_lowalpha(n, l, m, mu, M, a)
    end
    if alpha >= alpha_high
        return gamma_wkb(mu, M, a; kappa=kappa)
    end

    gamma_low = super_gamma_lowalpha(n, l, m, mu, M, a)
    gamma_hig = gamma_wkb(mu, M, a; kappa=kappa)
    w = (alpha - alpha_low) / (alpha_high - alpha_low)
    return exp((1.0 - w) * log(gamma_low) + w * log(gamma_hig))
end

function gamma_t_6g_to_5g(alpha, omega_tr, M)
    C = (2. ^28 * 3. ^4 * 5. ^5) / (11. ^22 * pi)
    C *= (32 * pi / 15.0)
    P_t = C * G * alpha^12 / r_g(M)^4
    return P_t / omega_tr
end

function dNedt(Ne, Ng, gamma_sr_yr, gamma_t_yr)
    return gamma_sr_yr * Ne - gamma_t_yr * Ne * Ng
end

function dNgdt(Ne, Ng, gamma_sr_yr, gamma_t_yr)
    return gamma_sr_yr * Ng + gamma_t_yr * Ne * Ng
end

function strain_envelope(dist_GeV_inv, alpha, ne, ng, mu, gamma_t_GeV, Ng, Ne)
    omega_tr = 0.5 * mu * alpha^2 * ((1.0 / ng^2) - (1.0 / ne^2))
    radicand = 4.0 * G / (dist_GeV_inv^2 * omega_tr) * gamma_t_GeV * Ng * Ne
    return sqrt(max(radicand, 0.0))
end

function linear_interp(x_new, x, y)
    result = similar(y, promote_type(eltype(y), Float64), length(x_new))
    n = length(x)

    if n == 0
        return result
    end

    if n == 1
        fill!(result, y[1])
        return result
    end

    for (idx, xi) in pairs(x_new)
        if xi <= x[1]
            result[idx] = y[1]
            continue
        end
        if xi >= x[end]
            result[idx] = y[end]
            continue
        end

        upper = searchsortedfirst(x, xi)
        lower = upper - 1
        x0 = x[lower]
        x1 = x[upper]
        y0 = y[lower]
        y1 = y[upper]
        t = (xi - x0) / (x1 - x0)
        result[idx] = y0 + (y1 - y0) * t
    end

    return result
end

function fftfreq(n::Int, d)
    result = Vector{Float64}(undef, n)
    half = div(n, 2)
    scale = 1.0 / (n * d)

    for i in 1:n
        if i <= half
            result[i] = (i - 1) * scale
        else
            result[i] = (i - 1 - n) * scale
        end
    end

    return result
end

function fft_continuous_from_time(time_yr, h_t; n_fft=2^20)
    time_yr = Float64.(collect(time_yr))
    h_t = ComplexF64.(collect(h_t))
    t_s = time_yr .* S_PER_YR
    t_uniform = range(t_s[1], t_s[end], length=n_fft)

    h_uniform = linear_interp(collect(t_uniform), t_s, real.(h_t))
    h_uniform = complex.(h_uniform, linear_interp(collect(t_uniform), t_s, imag.(h_t)))
    h_uniform .-= mean(h_uniform)

    dt = t_uniform[2] - t_uniform[1]
    H_f = fft(h_uniform) .* dt
    f_Hz = fftfreq(n_fft, dt)
    return f_Hz, H_f
end

function lorentzian(f, A, f0, gamma, C)
    return A ./ ((f .- f0).^2 .+ gamma^2) .+ C
end

function log_lorentzian(f, logA, f0, loggamma, logC)
    A = 10.0^logA
    gamma = 10.0^loggamma
    C = 10.0^logC
    return log10.(max.(A ./ ((f .- f0).^2 .+ gamma^2) .+ C, 1e-300))
end

function log_lorentzian(f, p::AbstractVector)
    return log_lorentzian(f, p[1], p[2], p[3], p[4])
end

function _numerical_jacobian(model, x, p, y0)
    n_params = length(p)
    n_points = length(y0)
    jacobian = zeros(eltype(y0), n_points, n_params)

    for j in 1:n_params
        step = max(abs(p[j]) * 1e-6, 1e-8)
        p_step = copy(p)
        p_step[j] += step
        y_step = model(x, p_step)
        jacobian[:, j] = (y_step .- y0) ./ step
    end

    return jacobian
end

function _least_squares_fit(model, x, y, p0; maxIter=200, tol=1e-10)
    p = Float64.(collect(p0))
    step_sizes = [
        max(abs(p[1]) * 0.25, 0.5),
        max(abs(p[2]) * 0.02, 1.0),
        max(abs(p[3]) * 0.25, 0.5),
        max(abs(p[4]) * 0.25, 0.5),
    ]

    best_error = sum(abs2, model(x, p) .- y)

    for _ in 1:maxIter
        improved = false

        for j in eachindex(p)
            for direction in (-1.0, 1.0)
                candidate = copy(p)
                candidate[j] += direction * step_sizes[j]
                candidate_error = sum(abs2, model(x, candidate) .- y)

                if isfinite(candidate_error) && candidate_error < best_error
                    p = candidate
                    best_error = candidate_error
                    improved = true
                end
            end
        end

        if !improved
            step_sizes .*= 0.5
            if maximum(step_sizes) < tol
                break
            end
        end
    end

    return p
end

function fit_lorentzian_to_fft(time_yr, h_t; n_fft=2^20, n_top=400)
    f_Hz, H_f = fft_continuous_from_time(time_yr, h_t; n_fft=n_fft)

    mask_pos = f_Hz .> 0
    f_pos = f_Hz[mask_pos]
    H_pos = H_f[mask_pos]
    abs_H_pos = abs.(H_pos)

    idx_sorted = sortperm(abs_H_pos, rev=true)
    n_top = min(n_top, length(idx_sorted))
    f_fit = f_pos[idx_sorted[1:n_top]]
    H_fit = abs_H_pos[idx_sorted[1:n_top]]

    sort_idx = sortperm(f_fit)
    f_fit = f_fit[sort_idx]
    H_fit = H_fit[sort_idx]

    f0_guess = f_fit[argmax(H_fit)]
    df = median(diff(sort(f_pos)))
    A0 = (maximum(H_fit) - minimum(H_fit)) * df^2
    gamma0 = df * 5.0
    C0 = minimum(H_fit)

    p0 = [log10(max(A0, 1e-40)), f0_guess, log10(max(gamma0, 1e-40)), log10(max(C0, 1e-40))]
    ydata = log10.(max.(H_fit, 1e-300))

    popt = _least_squares_fit(log_lorentzian, f_fit, ydata, p0; maxIter=20000)
    logA_fit, f0_fit, loggamma_fit, logC_fit = popt
    A_fit = 10.0^logA_fit
    gamma_fit = 10.0^loggamma_fit
    C_fit = 10.0^logC_fit

    L_full = A_fit ./ ((f_pos .- f0_fit).^2 .+ gamma_fit^2)

    @printf("Fitted centre frequency f0 = %.6e Hz\n", f0_fit)
    @printf("Fitted width gamma         = %.6e Hz\n", gamma_fit)
    @printf("Estimated FWHM ≃ 2·gamma   = %.6e Hz\n", 2 * gamma_fit)
    @printf("Fitted floor C             = %.3e\n", C_fit)

    return f_pos, H_pos, L_full, (A_fit, f0_fit, gamma_fit, C_fit)
end

end

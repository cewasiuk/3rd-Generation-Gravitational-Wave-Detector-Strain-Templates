if !isdefined(@__MODULE__, :ChrisWF)
    include("constants.jl")
end
using .ChrisWF
using Printf
using SpecialFunctions: expint

function _mp_E1_vec(z)
    return expint.(z)
end

function omega_ann(mua, alpha, n)
    return 2.0 * mua * (1.0 - alpha^2 / (2.0 * n^2))
end

function gamma_ann(l, alpha, M_GeV)
    p = l == 1 ? 17 : 4 * l + 1
    return G * 1e-10 / r_g(M_GeV)^3 * (((alpha / l) * 0.5)^p + ((alpha / l) * 0.5)^(p + 1))
end

function h_ann(delta_f_Hz, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)
    M_GeV = M_solar * M_sun
    r = r_kpc * kpc_to_GeV
    Gamma = gamma_ann(l, alpha, M_GeV)
    N_max = 10.0^76 * (M_solar / 10.0)^2

    delta_E = (collect(delta_f_Hz)) .* Hz_to_GeV
    z = im .* delta_E ./ (Gamma * N_max)

    pref = (1.0 / (2.0 * pi)) * sqrt(4.0 * G / (Gamma * r^2 * omega_a_GeV))

    Ei_part = _mp_E1_vec(z)
    Ei = exp.(im .* delta_E ./ (Gamma * N_max)) .* Ei_part

    return pref .* exp.(z) .* Ei .* h_GeV_s
end

function h_pcr_ann(f_grid_Hz, mua, M_solar, n, l, alpha, r_kpc, iota, phase)
    omega_a_GeV = omega_ann(mua, alpha, n)
    f_a_Hz = omega_a_GeV * GeV_to_Hz

    f_grid_Hz = collect(f_grid_Hz)
    delta_f_minus = f_grid_Hz .- f_a_Hz
    delta_f_plus = .-(f_grid_Hz .- f_a_Hz)

    hmin = h_ann(delta_f_minus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)
    hp = h_ann(delta_f_plus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)

    c = cos(iota)
    h_plus = (1.0 + c^2) / 4.0 * (exp(im * phase) * hmin + exp(-im * phase) * hp)
    h_cross = (c / (2.0im)) * (exp(im * phase) * hmin - exp(-im * phase) * hp)

    return h_plus, h_cross, f_grid_Hz, f_a_Hz
end

function h_pcr_ann_no_iota(f_grid_Hz, mua, M_solar, n, l, alpha, r_kpc, phase)
    omega_a_GeV = omega_ann(mua, alpha, n)
    f_a_Hz = omega_a_GeV * GeV_to_Hz

    f_grid_Hz = collect(f_grid_Hz)
    delta_f_minus = f_grid_Hz .- f_a_Hz
    delta_f_plus = .-(f_grid_Hz .- f_a_Hz)

    hmin = h_ann(delta_f_minus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)
    hp = h_ann(delta_f_plus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV)

    h_plus = 0.5 .* (exp(im * phase) * hmin + exp(-im * phase) * hp)
    h_cross = (1.0 / (2.0im)) .* (exp(im * phase) * hmin - exp(-im * phase) * hp)

    return h_plus, h_cross, f_grid_Hz, f_a_Hz
end

function amp_phase_pcr_ann_no_iota(f_grid_Hz, mua, M_solar, n, l, alpha, r_kpc, phase)
    omega_a_GeV = omega_ann(mua, alpha, n)
    f_a_Hz = omega_a_GeV * GeV_to_Hz

    f_grid_Hz = collect(f_grid_Hz)
    delta_f_minus = f_grid_Hz .- f_a_Hz

    h = exp(im * phase) .* h_ann(
        delta_f_minus, mua, M_solar, n, l, alpha, r_kpc, omega_a_GeV,
    )

    amp = abs.(h)
    phase_h = angle.(h)

    return amp, phase_h, f_grid_Hz, f_a_Hz
end

function iso_gatom_ann_strain(; 
    M_solar=3.1e-4,
    mua=2e-16,
    n=4,
    l=nothing,
    alpha=nothing,
    distance_kpc=1.0,
    iota=0.0,
    phase=0.0,
    f_min_Hz=1.0,
    f_max_Hz=1.0e12,
    n_f=50000,
    verbose=true,
)
    if l === nothing
        l = n - 1
    end

    if alpha === nothing
        alpha = G * (M_solar * M_sun) * mua
    end

    f_grid_Hz = 10.0 .^ range(log10(f_min_Hz), log10(f_max_Hz), length=n_f)

    h_plus, h_cross, f_grid_Hz, f_line_Hz = h_pcr_ann(
        f_grid_Hz, mua, M_solar, n, l, alpha, distance_kpc, iota, phase,
    )

    h_c = 2.0 .* f_grid_Hz .* sqrt.(abs.(h_plus).^2 .+ abs.(h_cross).^2)

    if verbose
        @printf("alpha           ≈ %.3e\n", alpha)
        @printf("r_g(M)          ≈ %.4e GeV^-1\n", r_g(M_solar * M_sun))
        @printf("Annihilation f_ann ≈ %.3e Hz\n", f_line_Hz)
        @printf("max |h_plus|    ≈ %.3e\n", maximum(abs.(h_plus)))
        @printf("max h_c         ≈ %.3e\n", maximum(h_c))
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

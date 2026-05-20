if !isdefined(Main, :ChrisWF)
    include("constants.jl")
end
using .ChrisWF
using Printf

function iso_gatom_level_tr_strain(; 
    M_solar=1e-6,
    a_spin=0.999999,
    alpha=1.0,
    ne=6,
    ng=5,
    m=nothing,
    distance_kpc=10.0,
    N_e0=1.0,
    N_g0=1.0,
    n_time=100000,
    n_fft=2^20,
    n_top=400,
    verbose=true,
)
    l = ng - 1
    if m === nothing
        m = l
    end

    M = M_solar * M_sun
    mu_a = alpha / (G * M)

    gamma_sre_GeV = super_gamma(ne, l, m, mu_a, M, a_spin)
    gamma_srg_GeV = super_gamma(ng, l, m, mu_a, M, a_spin)
    gamma_sre_yr = gamma_sre_GeV * GeV_to_yrinv
    gamma_srg_yr = gamma_srg_GeV * GeV_to_yrinv

    omega_tr = 0.5 * mu_a * alpha^2 * ((1.0 / ng^2) - (1.0 / ne^2))
    gamma_t_ne_GeV = gamma_t_6g_to_5g(alpha, omega_tr, M)
    gamma_t_ne_yr = gamma_t_ne_GeV * GeV_to_yrinv

    println()
    println(" alpha    = ", alpha)
    println(" omega_tr = ", omega_tr)
    @printf(" r_g(M)   = %.4e GeV^-1\n", r_g(M))
    @printf(" M        = %.4e GeV\n", M)

    dist_GeV_inv = distance_kpc * kpc_to_GeV

    dt_arr = 10.0 .^ range(-8, 1, length=n_time)
    yr = 0.0
    Ne = N_e0
    Ng = N_g0
    omega_sr = mu_a

    time_yr = Float64[]
    N_e_hist = Float64[]
    N_g_hist = Float64[]
    h_hist = Float64[]

    for j in dt_arr
        r_h = r_plus(M, a_spin)
        Omega_H = Omega_plus(a_spin, r_h)

        if m * Omega_H > omega_sr
            dNe = dNedt(Ne, Ng, gamma_sre_yr, gamma_t_ne_yr)
            dNg = dNgdt(Ne, Ng, gamma_srg_yr, gamma_t_ne_yr)

            Ne_next = Ne + dNe * j
            Ng_next = Ng + dNg * j

            if Ne_next < 0 || Ng_next < 0
                break
            end

            Ne = Ne_next
            Ng = Ng_next
            yr += j

            h_val = strain_envelope(
                dist_GeV_inv, alpha, ne, ng,
                mu_a, gamma_t_ne_GeV, Ng, Ne,
            )

            push!(N_e_hist, Ne)
            push!(N_g_hist, Ng)
            push!(h_hist, h_val)
            push!(time_yr, yr)
        else
            break
        end
    end

    time_yr = collect(time_yr)
    N_e_hist = collect(N_e_hist)
    N_g_hist = collect(N_g_hist)
    h_hist = Float64.(collect(h_hist))

    if verbose && !isempty(time_yr)
        @printf("Max N_e      ≈ %.3e\n", maximum(N_e_hist))
        @printf("Max N_g      ≈ %.3e\n", maximum(N_g_hist))
        @printf("Max h(t)     ≈ %.3e\n", maximum(abs.(h_hist)))
        @printf("omega_tr     ≈ %.3e GeV\n", omega_tr)
        @printf("gamma_sre_yr ≈ %.3e 1/yr\n", gamma_sre_yr)
        @printf("gamma_srg_yr ≈ %.3e 1/yr\n", gamma_srg_yr)
        @printf("gamma_t_yr   ≈ %.3e 1/yr\n", gamma_t_ne_yr)
        @printf("gamma_srg    ≈ %.3e 1/GeV\n", gamma_srg_GeV)
        @printf("Mu_a         ≈ %.3e GeV\n", mu_a)
    end

    f_pos, H_pos, L_full, lorentz_params = fit_lorentzian_to_fft(
        time_yr, h_hist; n_fft=n_fft, n_top=n_top,
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

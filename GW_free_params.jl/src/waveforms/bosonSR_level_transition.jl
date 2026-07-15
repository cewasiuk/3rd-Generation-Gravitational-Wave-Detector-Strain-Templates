module ChrisLevelTransition
include(joinpath(@__DIR__, "..", "..", "..", "chris", "julia_wf", "level_transition.jl"))
end

_numeric_value(x) = x
_numeric_value(x::ForwardDiff.Dual) = ForwardDiff.value(x)

const _BosonSRLevelModels = Union{BosonSR_level,BosonSR_level_FFT,BosonSR_level_num_der}

function _level_transition_frequency(model::_BosonSRLevelModels, M_solar, alpha)
    M = _numeric_value(M_solar) * ChrisLevelTransition.ChrisWF.M_sun
    alpha_val = _numeric_value(alpha)
    mu_a = alpha_val / (ChrisLevelTransition.ChrisWF.G * M)
    omega_tr = 0.5 * mu_a * alpha_val^2 * ((1.0 / model.ng^2) - (1.0 / model.ne^2))
    return omega_tr * ChrisLevelTransition.ChrisWF.GeV_to_Hz
end

function _level_lorentzian_params(model::_BosonSRLevelModels, M_solar, alpha, dL)
    dL_val = _numeric_value(dL)
    level_args = (;
        M_solar = _numeric_value(M_solar),
        a_spin = model.a_spin,
        alpha = _numeric_value(alpha),
        ne = model.ne,
        ng = model.ng,
        m = model.m,
        distance_kpc = dL_val * 1.0e6,
        N_e0 = model.N_e0,
        N_g0 = model.N_g0,
        n_time = model.n_time,
        n_fft = model.n_fft,
        n_top = model.n_top,
        verbose = model.verbose,
    )
    level = if model.verbose
        ChrisLevelTransition.iso_gatom_level_tr_strain(; level_args...)
    else
        redirect_stdout(devnull) do
            ChrisLevelTransition.iso_gatom_level_tr_strain(; level_args...)
        end
    end

    return level["lorentz_params"], level["f_pos"], level["L_full"], dL_val
end

function _level_spectrum(model::BosonSR_level, f::AbstractVector, M_solar, alpha, dL)
    _, f_pos, L_full, dL_val = _level_lorentzian_params(model, M_solar, alpha, dL)
    if isempty(f_pos)
        return zeros(eltype(f), length(f))
    end

    spectrum = linear_interpolation(f_pos, L_full, extrapolation_bc = 0.0)
    return spectrum(f) .* (dL_val / dL)
end

function _level_positive_sqrt(x)
    return x > zero(x) ? sqrt(x) : zero(x)
end

function _level_time_series(model::_BosonSRLevelModels, M_solar, alpha, dL)
    l = model.ng - 1
    m = model.m === nothing ? l : model.m

    M = M_solar * ChrisLevelTransition.ChrisWF.M_sun
    mu_a = alpha / (ChrisLevelTransition.ChrisWF.G * M)

    gamma_sre_GeV = ChrisLevelTransition.ChrisWF.super_gamma(model.ne, l, m, mu_a, M, model.a_spin)
    gamma_srg_GeV = ChrisLevelTransition.ChrisWF.super_gamma(model.ng, l, m, mu_a, M, model.a_spin)
    gamma_sre_yr = gamma_sre_GeV * ChrisLevelTransition.ChrisWF.GeV_to_yrinv
    gamma_srg_yr = gamma_srg_GeV * ChrisLevelTransition.ChrisWF.GeV_to_yrinv

    omega_tr = 0.5 * mu_a * alpha^2 * ((1.0 / model.ng^2) - (1.0 / model.ne^2))
    gamma_t_ne_GeV = ChrisLevelTransition.ChrisWF.gamma_t_6g_to_5g(alpha, omega_tr, M)
    gamma_t_ne_yr = gamma_t_ne_GeV * ChrisLevelTransition.ChrisWF.GeV_to_yrinv

    dist_GeV_inv = dL * 1.0e6 * ChrisLevelTransition.ChrisWF.kpc_to_GeV
    dt_arr = 10.0 .^ range(-8, 1, length = model.n_time)

    value_type = typeof(M + alpha + dL)
    seed = zero(M + alpha + dL)
    yr = seed
    Ne = seed + model.N_e0
    Ng = seed + model.N_g0
    omega_sr = mu_a

    time_yr = value_type[]
    h_hist = value_type[]

    for dt in dt_arr
        r_h = ChrisLevelTransition.ChrisWF.r_plus(M, model.a_spin)
        Omega_H = ChrisLevelTransition.ChrisWF.Omega_plus(model.a_spin, r_h)

        if m * Omega_H > omega_sr
            dNe = ChrisLevelTransition.ChrisWF.dNedt(Ne, Ng, gamma_sre_yr, gamma_t_ne_yr)
            dNg = ChrisLevelTransition.ChrisWF.dNgdt(Ne, Ng, gamma_srg_yr, gamma_t_ne_yr)

            Ne_next = Ne + dNe * dt
            Ng_next = Ng + dNg * dt

            if Ne_next < zero(Ne_next) || Ng_next < zero(Ng_next)
                break
            end

            Ne = Ne_next
            Ng = Ng_next
            yr += dt

            omega_tr_local = 0.5 * mu_a * alpha^2 * ((1.0 / model.ng^2) - (1.0 / model.ne^2))
            radicand = 4.0 * ChrisLevelTransition.ChrisWF.G / (dist_GeV_inv^2 * omega_tr_local) * gamma_t_ne_GeV * Ng * Ne
            h_val = _level_positive_sqrt(radicand)

            push!(time_yr, yr)
            push!(h_hist, h_val)
        else
            break
        end
    end

    return time_yr, h_hist
end

function _level_trapezoid_weights(t)
    n = length(t)
    if n == 0
        return eltype(t)[]
    elseif n == 1
        return [zero(eltype(t))]
    end

    weights = Vector{eltype(t)}(undef, n)
    weights[1] = (t[2] - t[1]) / 2
    weights[end] = (t[end] - t[end - 1]) / 2
    for idx in 2:(n - 1)
        weights[idx] = (t[idx + 1] - t[idx - 1]) / 2
    end
    return weights
end

function _level_spectrum(model::BosonSR_level_FFT, f::AbstractVector, M_solar, alpha, dL)
    time_yr, h_hist = _level_time_series(model, M_solar, alpha, dL)
    if length(time_yr) < 2
        sample = zero(M_solar + alpha + dL)
        return [complex(sample) for _ in f]
    end

    t_s = time_yr .* ChrisLevelTransition.ChrisWF.S_PER_YR
    weights = _level_trapezoid_weights(t_s)

    return [
        sum(@. weights * h_hist * exp(-2.0im * pi * fi * t_s))
        for fi in f
    ]
end

function _level_fit_params_for_num_der(model::BosonSR_level_num_der, M_solar, alpha, dL)
    params, _, _, _ = _level_lorentzian_params(model, M_solar, alpha, dL)
    return collect(params)
end

function _level_fit_param_jacobian(model::BosonSR_level_num_der, M0, alpha0, dL0)
    step_M = max(abs(M0) * model.fd_rel_step, model.fd_abs_step)
    step_alpha = max(abs(alpha0) * model.fd_rel_step, model.fd_abs_step)

    p_M_plus = _level_fit_params_for_num_der(model, M0 + step_M, alpha0, dL0)
    p_M_minus = _level_fit_params_for_num_der(model, M0 - step_M, alpha0, dL0)
    p_alpha_plus = _level_fit_params_for_num_der(model, M0, alpha0 + step_alpha, dL0)
    p_alpha_minus = _level_fit_params_for_num_der(model, M0, alpha0 - step_alpha, dL0)

    jac = Matrix{Float64}(undef, length(p_M_plus), 2)
    jac[:, 1] = (p_M_plus .- p_M_minus) ./ (2.0 * step_M)
    jac[:, 2] = (p_alpha_plus .- p_alpha_minus) ./ (2.0 * step_alpha)
    return jac
end

function _level_spectrum(model::BosonSR_level_num_der, f::AbstractVector, M_solar, alpha, dL)
    M0 = _numeric_value(M_solar)
    alpha0 = _numeric_value(alpha)
    dL0 = _numeric_value(dL)

    p0 = _level_fit_params_for_num_der(model, M0, alpha0, dL0)
    jac = _level_fit_param_jacobian(model, M0, alpha0, dL0)
    delta_M = M_solar - M0
    delta_alpha = alpha - alpha0
    p = p0 .+ jac[:, 1] .* delta_M .+ jac[:, 2] .* delta_alpha

    A, f0, gamma, _ = p
    return @. A / ((f - f0)^2 + gamma^2) * (dL0 / dL)
end

function Pol(model::_BosonSRLevelModels,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    M_solar, alpha = intrinsic_param
    amp = _level_spectrum(model, f, M_solar, alpha, dL)
    hp = @. 0.5 * (1.0 + cos(iota)^2) * amp
    hc = @. 1im * cos(iota) * amp
    return [hp, hc]
end

function Pol(model::_BosonSRLevelModels,
    f::AbstractVector,
    M_solar,
    alpha,
    dL,
    iota,
)
    return Pol(model, f, (M_solar, alpha), dL, iota)
end

function PolAbs(model::_BosonSRLevelModels,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    h_plus, h_cross = Pol(model, f, intrinsic_param, dL, iota)
    return [abs.(h_plus), abs.(h_cross)]
end

function PolAbs(model::_BosonSRLevelModels,
    f::AbstractVector,
    M_solar,
    alpha,
    dL,
    iota,
)
    return PolAbs(model, f, (M_solar, alpha), dL, iota)
end

function Phi(model::_BosonSRLevelModels,
    f::AbstractVector,
    intrinsic_param::Tuple,
    iota,
)
    return zeros(eltype(f), length(f))
end

function Phi(model::_BosonSRLevelModels,
    f::AbstractVector,
    M_solar,
    alpha,
    iota,
)
    return Phi(model, f, (M_solar, alpha), iota)
end

function Ampl(model::_BosonSRLevelModels,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    h_plus_abs, h_cross_abs = PolAbs(model, f, intrinsic_param, dL, iota)
    return sqrt.(h_plus_abs.^2 .+ h_cross_abs.^2)
end

function Ampl(model::_BosonSRLevelModels,
    f::AbstractVector,
    M_solar,
    alpha,
    dL,
    iota,
)
    return Ampl(model, f, (M_solar, alpha), dL, iota)
end

function _fcut(model::_BosonSRLevelModels, intrinsic_param::Tuple)
    return _fcut(model, intrinsic_param...)
end

function _fcut(model::_BosonSRLevelModels, M_solar, alpha)
    return max(2.0, 1.1 * _level_transition_frequency(model, M_solar, alpha))
end

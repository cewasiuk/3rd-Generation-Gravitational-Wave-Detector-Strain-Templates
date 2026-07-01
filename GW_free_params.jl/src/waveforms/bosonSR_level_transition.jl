module ChrisLevelTransition
include(joinpath(@__DIR__, "..", "..", "..", "chris", "julia_wf", "level_transition.jl"))
end

_numeric_value(x) = x
_numeric_value(x::ForwardDiff.Dual) = ForwardDiff.value(x)

function _level_transition_frequency(model::BosonSR_level, M_solar, alpha)
    M = _numeric_value(M_solar) * ChrisLevelTransition.ChrisWF.M_sun
    alpha_val = _numeric_value(alpha)
    mu_a = alpha_val / (ChrisLevelTransition.ChrisWF.G * M)
    omega_tr = 0.5 * mu_a * alpha_val^2 * ((1.0 / model.ng^2) - (1.0 / model.ne^2))
    return omega_tr * ChrisLevelTransition.ChrisWF.GeV_to_Hz
end

function _level_spectrum(model::BosonSR_level, f::AbstractVector, M_solar, alpha, dL)
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

    f_pos = level["f_pos"]
    L_full = level["L_full"]
    if isempty(f_pos)
        return zeros(eltype(f), length(f))
    end

    spectrum = linear_interpolation(f_pos, L_full, extrapolation_bc = 0.0)
    return spectrum(f) .* (dL_val / dL)
end

function Pol(model::BosonSR_level,
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

function Pol(model::BosonSR_level,
    f::AbstractVector,
    M_solar,
    alpha,
    dL,
    iota,
)
    return Pol(model, f, (M_solar, alpha), dL, iota)
end

function PolAbs(model::BosonSR_level,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    h_plus, h_cross = Pol(model, f, intrinsic_param, dL, iota)
    return [abs.(h_plus), abs.(h_cross)]
end

function PolAbs(model::BosonSR_level,
    f::AbstractVector,
    M_solar,
    alpha,
    dL,
    iota,
)
    return PolAbs(model, f, (M_solar, alpha), dL, iota)
end

function Phi(model::BosonSR_level,
    f::AbstractVector,
    intrinsic_param::Tuple,
    iota,
)
    return zeros(eltype(f), length(f))
end

function Phi(model::BosonSR_level,
    f::AbstractVector,
    M_solar,
    alpha,
    iota,
)
    return Phi(model, f, (M_solar, alpha), iota)
end

function Ampl(model::BosonSR_level,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    h_plus_abs, h_cross_abs = PolAbs(model, f, intrinsic_param, dL, iota)
    return sqrt.(h_plus_abs.^2 .+ h_cross_abs.^2)
end

function Ampl(model::BosonSR_level,
    f::AbstractVector,
    M_solar,
    alpha,
    dL,
    iota,
)
    return Ampl(model, f, (M_solar, alpha), dL, iota)
end

function _fcut(model::BosonSR_level, intrinsic_param::Tuple)
    return _fcut(model, intrinsic_param...)
end

function _fcut(model::BosonSR_level, M_solar, alpha)
    return max(2.0, 1.1 * _level_transition_frequency(model, M_solar, alpha))
end

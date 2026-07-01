module HenryBinaryTransition
include(joinpath(@__DIR__, "..", "..", "..", "henry", "julia_wf", "functions.jl"))
end

const _HENRY_HZ_PER_EV = 2.417987242e14
const _HENRY_GPC_TO_EV_INV = 1.5637381227081487e32

function _henry_frequency_eV(f::AbstractVector)
    return f ./ _HENRY_HZ_PER_EV
end

function _henry_distance_eV_inv(dL)
    return _numeric_value(dL) * _HENRY_GPC_TO_EV_INV
end

function _henry_central_frequency_Hz(model::BosonSR_binary, q, M_solar, alpha)
    M = _numeric_value(M_solar) * HenryBinaryTransition.HenryWF.Msol_in_eV
    Omega0 = HenryBinaryTransition.HenryWF.Omega0_binary_natural_unit(
        model.m_i,
        _numeric_value(alpha),
        HenryBinaryTransition.HenryWF.G,
        M,
        model.n,
        model.l_i,
    )
    return HenryBinaryTransition.HenryWF.fc_from_Omega0(Omega0) * _HENRY_HZ_PER_EV
end

function _henry_plus_amplitude(model::BosonSR_binary,
    f::AbstractVector,
    q,
    M_solar,
    alpha,
    Gamma_abs,
    dL,
)
    M_solar_val = _numeric_value(M_solar)
    alpha_val = _numeric_value(alpha)
    q_val = _numeric_value(q)
    dL_val = _numeric_value(dL)
    M = M_solar_val * HenryBinaryTransition.HenryWF.Msol_in_eV
    Omega0 = HenryBinaryTransition.HenryWF.Omega0_binary_natural_unit(
        model.m_i,
        alpha_val,
        HenryBinaryTransition.HenryWF.G,
        M,
        model.n,
        model.l_i,
    )
    eta = HenryBinaryTransition.HenryWF.eta_parameter(alpha_val, q_val, M_solar_val)
    h_of_iota = HenryBinaryTransition.HenryWF.htilde_plus(
        _henry_frequency_eV(f),
        M,
        _henry_distance_eV_inv(dL_val),
        alpha_val,
        Omega0,
        q_val,
        model.m_i,
        model.m_f,
        eta,
        _numeric_value(Gamma_abs);
        use_z_scaling = model.use_z_scaling,
        numerical_qc = model.numerical_qc,
    )
    return h_of_iota(0.0) .* (dL_val / dL)
end

function Pol(model::BosonSR_binary,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    q, M_solar, alpha, Gamma_abs = intrinsic_param
    h0 = _henry_plus_amplitude(model, f, q, M_solar, alpha, Gamma_abs, dL)
    h_plus = @. 0.5 * (1.0 + cos(iota)^2) * h0
    h_cross = zeros(Complex{eltype(h_plus)}, length(h_plus))
    return [h_plus, h_cross]
end

function Pol(model::BosonSR_binary,
    f::AbstractVector,
    q,
    M_solar,
    alpha,
    Gamma_abs,
    dL,
    iota,
)
    return Pol(model, f, (q, M_solar, alpha, Gamma_abs), dL, iota)
end

function PolAbs(model::BosonSR_binary,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    h_plus, h_cross = Pol(model, f, intrinsic_param, dL, iota)
    return [abs.(h_plus), abs.(h_cross)]
end

function PolAbs(model::BosonSR_binary,
    f::AbstractVector,
    q,
    M_solar,
    alpha,
    Gamma_abs,
    dL,
    iota,
)
    return PolAbs(model, f, (q, M_solar, alpha, Gamma_abs), dL, iota)
end

function Phi(model::BosonSR_binary,
    f::AbstractVector,
    intrinsic_param::Tuple,
    iota,
)
    return zeros(eltype(f), length(f))
end

function Phi(model::BosonSR_binary,
    f::AbstractVector,
    q,
    M_solar,
    alpha,
    Gamma_abs,
    iota,
)
    return Phi(model, f, (q, M_solar, alpha, Gamma_abs), iota)
end

function Ampl(model::BosonSR_binary,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    h_plus_abs, h_cross_abs = PolAbs(model, f, intrinsic_param, dL, iota)
    return sqrt.(h_plus_abs.^2 .+ h_cross_abs.^2)
end

function Ampl(model::BosonSR_binary,
    f::AbstractVector,
    q,
    M_solar,
    alpha,
    Gamma_abs,
    dL,
    iota,
)
    return Ampl(model, f, (q, M_solar, alpha, Gamma_abs), dL, iota)
end

function _fcut(model::BosonSR_binary, intrinsic_param::Tuple)
    return _fcut(model, intrinsic_param...)
end

function _fcut(model::BosonSR_binary, q, M_solar, alpha, Gamma_abs)
    return max(2.0, 1.5 * _henry_central_frequency_Hz(model, q, M_solar, alpha))
end

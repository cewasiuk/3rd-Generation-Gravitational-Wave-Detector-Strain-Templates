module HenryBinaryTransition
include(joinpath(@__DIR__, "..", "..", "..", "henry", "julia_wf", "functions.jl"))
end

const _HENRY_HZ_PER_EV = 2.417987242e14
const _HENRY_GPC_TO_EV_INV = 1.5637381227081487e32

function _henry_frequency_eV(f::AbstractVector)
    return f ./ _HENRY_HZ_PER_EV
end

function _henry_distance_eV_inv(dL)
    return dL * _HENRY_GPC_TO_EV_INV
end

function _henry_use_numerical_qc(model::BosonSR_binary)
    return model.numerical_qc
end

function _henry_use_numerical_qc(model::BosonSR_binary_num)
    return true
end

function _henry_central_frequency_Hz(model::_BosonSRBinaryModels, q, M_solar, alpha)
    M = M_solar * HenryBinaryTransition.HenryWF.Msol_in_eV
    Omega0 = HenryBinaryTransition.HenryWF.Omega0_binary_natural_unit(
        model.m_i,
        alpha,
        HenryBinaryTransition.HenryWF.G,
        M,
        model.n,
        model.l_i,
    )
    return HenryBinaryTransition.HenryWF.fc_from_Omega0(Omega0) * _HENRY_HZ_PER_EV
end

function _henry_polarizations(model::_BosonSRBinaryModels,
    f::AbstractVector,
    q,
    M_solar,
    alpha,
    Gamma_abs,
    dL,
    iota,
)
    M = M_solar * HenryBinaryTransition.HenryWF.Msol_in_eV
    Omega0 = HenryBinaryTransition.HenryWF.Omega0_binary_natural_unit(
        model.m_i,
        alpha,
        HenryBinaryTransition.HenryWF.G,
        M,
        model.n,
        model.l_i,
    )
    eta = HenryBinaryTransition.HenryWF.eta_parameter(alpha, q, M_solar)
    h_plus_of_iota = HenryBinaryTransition.HenryWF.htilde_plus(
        _henry_frequency_eV(f),
        M,
        _henry_distance_eV_inv(dL),
        alpha,
        Omega0,
        q,
        model.m_i,
        model.m_f,
        eta,
        Gamma_abs;
        use_z_scaling = model.use_z_scaling,
        numerical_qc = _henry_use_numerical_qc(model),
    )
    h_cross_of_iota = HenryBinaryTransition.HenryWF.htilde_cross(
        _henry_frequency_eV(f),
        M,
        _henry_distance_eV_inv(dL),
        alpha,
        Omega0,
        q,
        model.m_i,
        model.m_f,
        eta,
        Gamma_abs;
        use_z_scaling = model.use_z_scaling,
        numerical_qc = _henry_use_numerical_qc(model),
    )
    return h_plus_of_iota(iota), h_cross_of_iota(iota)
end

function _henry_plus_amplitude(model::_BosonSRBinaryModels,
    f::AbstractVector,
    q,
    M_solar,
    alpha,
    Gamma_abs,
    dL,
)
    h_plus, _ = _henry_polarizations(model, f, q, M_solar, alpha, Gamma_abs, dL, 0.0)
    return h_plus
end

function Pol(model::_BosonSRBinaryModels,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    q, M_solar, alpha, Gamma_abs = intrinsic_param
    h_plus, h_cross = _henry_polarizations(model, f, q, M_solar, alpha, Gamma_abs, dL, iota)
    return [h_plus, h_cross]
end

function Pol(model::_BosonSRBinaryModels,
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

function PolAbs(model::_BosonSRBinaryModels,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    h_plus, h_cross = Pol(model, f, intrinsic_param, dL, iota)
    return [abs.(h_plus), abs.(h_cross)]
end

function PolAbs(model::_BosonSRBinaryModels,
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

function Phi(model::_BosonSRBinaryModels,
    f::AbstractVector,
    intrinsic_param::Tuple,
    iota,
)
    return zeros(eltype(f), length(f))
end

function Phi(model::_BosonSRBinaryModels,
    f::AbstractVector,
    q,
    M_solar,
    alpha,
    Gamma_abs,
    iota,
)
    return Phi(model, f, (q, M_solar, alpha, Gamma_abs), iota)
end

function Ampl(model::_BosonSRBinaryModels,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
)
    h_plus_abs, h_cross_abs = PolAbs(model, f, intrinsic_param, dL, iota)
    return sqrt.(h_plus_abs.^2 .+ h_cross_abs.^2)
end

function Ampl(model::_BosonSRBinaryModels,
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

function _fcut(model::_BosonSRBinaryModels, intrinsic_param::Tuple)
    return _fcut(model, intrinsic_param...)
end

function _fcut(model::_BosonSRBinaryModels, q, M_solar, alpha, Gamma_abs)
    return max(2.0, 1.5 * _henry_central_frequency_Hz(model, q, M_solar, alpha))
end

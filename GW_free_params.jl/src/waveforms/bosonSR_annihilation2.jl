module ChrisAnnihilation2
include(joinpath(@__DIR__, "..", "..", "..", "chris", "julia_wf", "annihilation2.jl"))
end

const _GPC_TO_AU = 1.0e6 * 2.0626480624709636e8

function _ann2_alpha(M_solar, mua)
    return ChrisAnnihilation2.alpha_of(M_solar, mua)
end

function _ann2_line_frequency(M_solar, mua; n=2)
    alpha = _ann2_alpha(M_solar, mua)
    return ChrisAnnihilation2.f_GW_Hz(mua, alpha, n)
end

function Pol(model::BosonSR_ann2,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
    phase=0.0;
    n=2,
    cloud_efficiency=ChrisAnnihilation2.CLOUD_EFFICIENCY_DEFAULT,
    M_a_GeV=nothing,
    a_star=0.0,
    T_obs_s=4.0 * 86400.0,
)
    M_solar, mua = intrinsic_param
    result = ChrisAnnihilation2.iso_gatom_ann_strain(;
        M_solar=M_solar,
        mua=mua,
        n=n,
        distance_au=dL * _GPC_TO_AU,
        iota=iota,
        phase=phase,
        cloud_efficiency=cloud_efficiency,
        M_a_GeV=M_a_GeV,
        a_star=a_star,
        T_obs_s=T_obs_s,
        verbose=false,
    )
    h_plus, h_cross = ChrisAnnihilation2.sampled_strain(
        f,
        result["f_line_Hz"],
        result["h_plus"],
        result["h_cross"],
        T_obs_s,
        result["tau_s"],
    )
    return [h_plus, h_cross]
end

function Pol(model::BosonSR_ann2,
    f::AbstractVector,
    M_solar,
    mua,
    dL,
    iota,
    phase=0.0;
    kwargs...
)
    return Pol(model, f, (M_solar, mua), dL, iota, phase; kwargs...)
end

function PolAbs(model::BosonSR_ann2,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
    phase=0.0;
    kwargs...
)
    h_plus, h_cross = Pol(model, f, intrinsic_param, dL, iota, phase; kwargs...)
    return [abs.(h_plus), abs.(h_cross)]
end

function PolAbs(model::BosonSR_ann2,
    f::AbstractVector,
    M_solar,
    mua,
    dL,
    iota,
    phase=0.0;
    kwargs...
)
    return PolAbs(model, f, (M_solar, mua), dL, iota, phase; kwargs...)
end

function Phi(model::BosonSR_ann2,
    f::AbstractVector,
    intrinsic_param::Tuple,
    iota,
    phase=0.0;
    kwargs...
)
    return zeros(eltype(f), length(f))
end

function Phi(model::BosonSR_ann2,
    f::AbstractVector,
    M_solar,
    mua,
    iota,
    phase=0.0;
    kwargs...
)
    return Phi(model, f, (M_solar, mua), iota, phase; kwargs...)
end

function Ampl(model::BosonSR_ann2,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
    phase=0.0;
    kwargs...
)
    h_plus_abs, h_cross_abs = PolAbs(model, f, intrinsic_param, dL, iota, phase; kwargs...)
    return sqrt.(h_plus_abs.^2 .+ h_cross_abs.^2)
end

function Ampl(model::BosonSR_ann2,
    f::AbstractVector,
    M_solar,
    mua,
    dL,
    iota,
    phase=0.0;
    kwargs...
)
    return Ampl(model, f, (M_solar, mua), dL, iota, phase; kwargs...)
end

function _fcut(model::BosonSR_ann2, intrinsic_param::Tuple)
    return _fcut(model, intrinsic_param...)
end

function _fcut(model::BosonSR_ann2, M_solar, mua)
    return max(2.0, 1.1 * _ann2_line_frequency(M_solar, mua))
end

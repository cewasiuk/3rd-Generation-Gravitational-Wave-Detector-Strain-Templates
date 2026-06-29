module ChrisAnnihilation
include(joinpath(@__DIR__, "..", "..", "..", "chris", "julia_wf", "annihilation.jl"))
end

function _ann_alpha(M_solar, mua)
    return ChrisAnnihilation.ChrisWF.G * (M_solar * ChrisAnnihilation.ChrisWF.M_sun) * mua
end

function _ann_line_frequency(M_solar, mua; n=4)
    alpha = _ann_alpha(M_solar, mua)
    omega_a_GeV = ChrisAnnihilation.omega_ann(mua, alpha, n)
    return omega_a_GeV * ChrisAnnihilation.ChrisWF.GeV_to_Hz
end

function Pol(model::BosonSR_ann,
    f::AbstractVector,
    intrinsic_param::Tuple,
    dL,
    iota,
    phase=0.0;
    n=4,
    l=n-1,
)
    M_solar, mua = intrinsic_param
    alpha = _ann_alpha(M_solar, mua)
    distance_kpc = dL * 1e6

    h_plus, h_cross, _, _ = ChrisAnnihilation.h_pcr_ann(
        f, mua, M_solar, n, l, alpha, distance_kpc, iota, phase,
    )

    return [h_plus, h_cross]
end

function Pol(model::BosonSR_ann,
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

function PolAbs(model::BosonSR_ann,
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

function PolAbs(model::BosonSR_ann,
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

function Phi(model::BosonSR_ann,
    f::AbstractVector,
    intrinsic_param::Tuple,
    iota,
    phase=0.0;
    kwargs...
)
    return zeros(eltype(f), length(f))
end

function Phi(model::BosonSR_ann,
    f::AbstractVector,
    M_solar,
    mua,
    iota,
    phase=0.0;
    kwargs...
)
    return Phi(model, f, (M_solar, mua), iota, phase; kwargs...)
end

function Ampl(model::BosonSR_ann,
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

function Ampl(model::BosonSR_ann,
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

function _fcut(model::BosonSR_ann, intrinsic_param::Tuple)
    return _fcut(model, intrinsic_param...)
end

function _fcut(model::BosonSR_ann,
    M_solar,
    mua,
)
    return max(2.0, 1.1 * _ann_line_frequency(M_solar, mua))
end

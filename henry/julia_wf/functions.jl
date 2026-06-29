module HenryWF

using DelimitedFiles
using LinearAlgebra
using Printf

# =======================================
# Natural Units
# =======================================
const Mp = 1.220890e28 # eV
const G = 1.0 / Mp^2

# =======================================
# Unit Conversion from BHSR
# =======================================
const eV_per_kg = 1.782622e-36
const mP_in_GeV = 1.220890e19 # GeV
const mP_in_eV = 1e9 * mP_in_GeV # eV

# Astro constants
const Msol_in_kg = 1.99841e30 # kg
const Msol_in_eV = Msol_in_kg / eV_per_kg # eV

const mP2_in_eVMsol = mP_in_eV * (mP_in_eV / Msol_in_eV) # eV Msol
const GNewton = 1.0 / mP2_in_eVMsol # eV^-1 Msol^-1

const _LANCZOS_G = 7.0
const _LANCZOS_COEFFS = (
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
)

export Mp, G, eV_per_kg, mP_in_GeV, mP_in_eV, Msol_in_kg, Msol_in_eV,
    mP2_in_eVMsol, GNewton,
    Omega0_binary_natural_unit, gamma_rates_natural_unit, eta_parameter,
    set_units, electron_mass, fac, Cslm, s_lambda_lm, sYlm,
    rg, alpha, r_plus, omegaHyperfine, omega0_bxzh, omega1_bxzh,
    gam_pq_bxzh, omega_nlm_bxzh, h_seidel, flm_seidel_2, flm_seidel_4,
    alm_approx, swsphericalh_A, M_matrix_elem, l_min, ells, M_matrix,
    sep_consts, angular_ev, cfunctions, alpha_n, beta_n, gamma_n,
    continued_fraction, root_equation, find_cf_root, cfm_bhsr_rates,
    a_tilde_crit, q_c, h0_from_params, gamma_rate, z_parameter,
    z_scaling_211_to_21m1, fc_from_Omega0, psi_plus, htilde_plus,
    compute_cloud_mass_numerical

function _gamma_lanczos(z::Complex)
    if real(z) < 0.5
        return pi / (sin(pi * z) * _gamma_lanczos(1.0 - z))
    end

    z1 = z - 1.0
    x = complex(_LANCZOS_COEFFS[1])
    for i in 2:length(_LANCZOS_COEFFS)
        x += _LANCZOS_COEFFS[i] / (z1 + (i - 1))
    end
    t = z1 + _LANCZOS_G + 0.5
    return sqrt(2.0 * pi) * t^(z1 + 0.5) * exp(-t) * x
end

_gamma_lanczos(x::Real) = _gamma_lanczos(complex(float(x)))

function _simpson(f, a, b)
    c = 0.5 * (a + b)
    return (b - a) * (f(a) + 4.0 * f(c) + f(b)) / 6.0
end

function _adaptive_simpson(f, a, b, eps, whole, depth)
    c = 0.5 * (a + b)
    left = _simpson(f, a, c)
    right = _simpson(f, c, b)
    delta = left + right - whole
    if depth <= 0 || abs(delta) <= 15.0 * eps
        return left + right + delta / 15.0
    end
    return _adaptive_simpson(f, a, c, eps / 2.0, left, depth - 1) +
           _adaptive_simpson(f, c, b, eps / 2.0, right, depth - 1)
end

function quad_integral(f, a, b; rtol=1e-8, maxdepth=22)
    if a == b
        return 0.0
    end
    whole = _simpson(f, a, b)
    return _adaptive_simpson(f, a, b, rtol, whole, maxdepth)
end

function quad_to_infinity(f, a; rtol=1e-8, maxdepth=24)
    eps = 1e-12
    g(t) = begin
        r = a + t / (1.0 - t)
        f(r) / (1.0 - t)^2
    end
    return quad_integral(g, 0.0, 1.0 - eps; rtol=rtol, maxdepth=maxdepth)
end

function _factorial_float(n::Integer)
    if n < 0
        throw(DomainError(n, "factorial is undefined for negative integers"))
    end
    return Float64(factorial(big(n)))
end

function _assoc_legendre(l::Integer, m::Integer, x)
    if m < 0 || m > l
        throw(ArgumentError("associated Legendre implementation expects 0 <= m <= l"))
    end
    pmm = 1.0
    if m > 0
        somx2 = sqrt(max(0.0, (1.0 - x) * (1.0 + x)))
        fact = 1.0
        for _ in 1:m
            pmm *= -fact * somx2
            fact += 2.0
        end
    end
    if l == m
        return pmm
    end
    pmmp1 = x * (2.0 * m + 1.0) * pmm
    if l == m + 1
        return pmmp1
    end
    pll = 0.0
    for ll in (m + 2):l
        pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m)
        pmm = pmmp1
        pmmp1 = pll
    end
    return pll
end

function spherical_harmonic(m::Integer, l::Integer, phi, theta)
    if abs(m) > l
        return 0.0 + 0.0im
    end
    if m < 0
        return (-1)^m * conj(spherical_harmonic(-m, l, phi, theta))
    end
    norm = sqrt((2.0 * l + 1.0) / (4.0 * pi) * _factorial_float(l - m) / _factorial_float(l + m))
    return norm * _assoc_legendre(l, m, cos(theta)) * exp(im * m * phi)
end

# =======================================
# Numerical Calculation of the Mixing Amplitude Eta
# =======================================

function Omega0_binary_natural_unit(m_i, alpha, G_in, M, n, l_i)
    return (64.0 * m_i * alpha^7) /
           (G_in * M * n^3 * 2.0 * l_i * (2.0 * l_i + 1.0) * (2.0 * l_i + 2.0) * (m_i^2 + 4.0 * alpha^2))
end

function gamma_rates_natural_unit(q, G_in, M, Omega0)
    return (96.0 / 5.0) * (q / (1.0 + q)^(1.0 / 3.0)) * (G_in * M * Omega0)^(5.0 / 3.0) * Omega0^2
end

function eta_parameter(alpha_in, q, mass)
    M = mass * Msol_in_eV

    n = 2
    l_star, m_star = 2, -2
    l_i, m_i = 1, 1
    l_f, m_f = 1, -1
    Delta_m = abs(m_f - m_i)

    Omega0_binary_natural = Omega0_binary_natural_unit(m_i, alpha_in, G, M, n, l_i)
    gamma_rates_natural_unit(q, G, M, Omega0_binary_natural)

    function angular_integrand(theta, phi)
        Y_star = spherical_harmonic(m_star, l_star, phi, theta)
        Y_i = spherical_harmonic(m_i, l_i, phi, theta)
        Y_f_conj = conj(spherical_harmonic(m_f, l_f, phi, theta))
        return real(Y_star * Y_i * Y_f_conj * sin(theta))
    end

    I_A = quad_integral(phi -> quad_integral(theta -> angular_integrand(theta, phi), 0.0, pi; rtol=1e-8), 0.0, 2.0 * pi; rtol=1e-8)

    r_c = G * M / alpha_in^2
    R_star = (G * M / Omega0_binary_natural^2)^(1.0 / 3.0)
    r_bar_max = R_star / r_c

    R_21(r_bar) = (1.0 / sqrt(24.0)) * r_bar * exp(-r_bar / 2.0)

    I_in = quad_integral(r_bar -> R_21(r_bar)^2 * r_bar^4, 0.0, r_bar_max; rtol=1e-8)
    I_out = quad_to_infinity(r_bar -> R_21(r_bar)^2 / r_bar, r_bar_max; rtol=1e-8)

    prefactor = sqrt(3.0 * pi / 10.0) * I_A
    term1 = (q * G * M * Omega0_binary_natural) / (alpha_in^3 * (1.0 + q)) * I_in
    term2 = (alpha_in^7 * q * (1.0 + q)^(2.0 / 3.0)) / (G * M * Omega0_binary_natural)^(7.0 / 3.0) * I_out
    eta = Omega0_binary_natural * abs(prefactor * (term1 + term2))
    return eta
end

# =======================================
# Numerical Calculation of the Cloud Mass q_c
# =======================================

function set_units(units, Mbh0)
    if units == "physical" || units == "physical+alpha"
        tunit = 4.920551932748678e-06
        Punit = 3.6283745e52
        dunit = 4.78691895e-20
        hbar = 1.19727031e-76
        muunit = units == "physical+alpha" ? 1.0 / Mbh0 : 7.48548859e9
    elseif units == "natural" || units == "natural+alpha"
        tunit = 1.0
        Punit = 1.0
        dunit = 1.0
        hbar = 1.0
        muunit = units == "natural+alpha" ? 1.0 / Mbh0 : 1.0
    else
        throw(ArgumentError("Invalid boson cloud waveform units"))
    end
    return (tunit, Punit, dunit, hbar, muunit)
end

function electron_mass(units)
    if units == "physical" || units == "physical+alpha"
        return 3.825076814891968e15
    elseif units == "natural" || units == "natural+alpha"
        return 4.18543e-23
    else
        throw(ArgumentError("Invalid boson cloud waveform units"))
    end
end

function fac(n::Integer)
    result = big(1)
    for i in 2:n
        result *= i
    end
    return result
end

function Cslm(s, l, m)
    return sqrt(l * l * (4.0 * l * l - 1.0) / ((l * l - m * m) * (l * l - s * s)))
end

function s_lambda_lm(s, l, m, x)
    Pm = (-0.5)^m
    if m != s
        Pm = Pm .* (1.0 .+ x).^((m - s) * 1.0 / 2.0)
    end
    if m != -s
        Pm = Pm .* (1.0 .- x).^((m + s) * 1.0 / 2.0)
    end
    Pm = Pm .* sqrt(Float64(fac(2 * m + 1)) / (4.0 * pi * Float64(fac(m + s)) * Float64(fac(m - s))))
    if l == m
        return Pm .* ones(size(x))
    end
    Pm1 = (x .+ s * 1.0 / (m + 1)) .* Cslm(s, m + 1, m) .* Pm
    if l == m + 1
        return Pm1
    else
        Pn = Pm1
        for n in (m + 2):l
            Pn = (x .+ s * m * 1.0 / (n * (n - 1.0))) .* Cslm(s, n, m) .* Pm1 .-
                 Cslm(s, n, m) * 1.0 / Cslm(s, n - 1, m) .* Pm
            Pm = Pm1
            Pm1 = Pn
        end
        return Pn
    end
end

function sYlm(ss, ll, mm, theta)
    l = ll
    m = mm
    s = ss
    if l < 0
        return 0
    end
    if abs(m) > l || l < abs(s)
        return 0
    end
    Pm = 1.0
    if abs(mm) < abs(ss)
        s = mm
        m = ss
        if isodd(m + s)
            Pm = -Pm
        end
    end
    if m < 0
        s = -s
        m = -m
        if isodd(m + s)
            Pm = -Pm
        end
    end
    return Pm .* s_lambda_lm(s, l, m, cos(theta))
end

function compute_cloud_mass_numerical(alpha_in, M, boson_mass; spin=0.99)
    throw(ArgumentError(
        "compute_cloud_mass_numerical is table/interpolator backed in the Python source. " *
        "Use htilde_plus(...; numerical_qc=false) for the analytical q_c path, or port " *
        "RelScalar/MatchedWaveform with native .npz interpolation before enabling this."
    ))
end

# =======================================
# Functions for calculating Relativitistic Gamma Rate from BHSR starts here
# =======================================

function rg(mbh::Real)
    return GNewton * mbh
end

function alpha(mu, mbh)
    return rg(mbh) .* mu
end

function r_plus(mbh::Real, astar::Real)
    return rg(mbh) * (1.0 + sqrt(1.0 - astar * astar))
end

function omegaHyperfine(mu::Real, mbh::Real, astar::Real, n::Integer, l::Integer, m::Integer)
    x = alpha(mu, mbh) / n
    x2 = x * x
    x4 = x2 * x2
    fine = 1.875 - 6.0 * n / (2 * l + 1)
    hyperfine = 8.0 * m * n * n * astar / (l * (2 * l + 1) * (2 * l + 2))
    return mu * (1.0 - 0.5 * x2 + fine * x4 + hyperfine * x * x4)
end

function omega0_bxzh(mu::Real, mbh::Real, n::Integer)
    n2 = n * n
    al = alpha(mu, mbh)
    al2 = al * al
    x = 2.0 * al2 / (n2 + 4.0 * al2 + n * sqrt(n2 + 8.0 * al2))
    return mu * sqrt(1.0 - x)
end

function omega1_bxzh(mu::Real, mbh::Real, n::Integer)
    om0 = omega0_bxzh(mu, mbh, n)
    if om0 > 0
        om02 = om0 * om0
        mu2 = mu * mu
        al = alpha(mu, mbh)
        al2 = al * al
        x = 1.0 + 4.0 * al2 * (2.0 * om02 / mu2 - 1.0) / (n * n)
        return (mu2 - om02) / (n * om0 * x)
    end
    return 0.0
end

function gam_pq_bxzh(p::Real, q, eps::Real, n::Integer, l::Integer)
    lp = l + eps
    ip = p * im
    twolp = 2.0 * lp
    g1 = _gamma_lanczos(twolp + 1.0)
    g2 = _gamma_lanczos(twolp + 2.0)
    g2n = _gamma_lanczos(twolp + 1.0 + n - l)
    gpmeps = _gamma_lanczos(1.0 + 2.0 * eps) * _gamma_lanczos(1.0 - 2.0 * eps)
    x1 = lp + 1.0 + ip
    x2 = sqrt(complex(q - p * p))
    gabs = abs(_gamma_lanczos(x1 + x2) * _gamma_lanczos(x1 - x2))
    gmix = _gamma_lanczos(1.0 - ip - eps + x2) * _gamma_lanczos(1.0 - ip - eps - x2)
    gmix *= _gamma_lanczos(1.0 + ip + eps + x2) * _gamma_lanczos(1.0 + ip + eps - x2)
    num = g2n * gpmeps * gabs * gabs * 2.0^(4.0 * lp + 2.0)
    denom = Float64(fac(n - l - 1)) * g1 * g1 * g2 * g2 * gmix
    return num / denom
end

function omega_nlm_bxzh(mu::Real, mbh::Real, astar::Real, n::Integer=2, l::Integer=1, m::Integer=1)
    om0 = omega0_bxzh(mu, mbh, n)
    om1 = omega1_bxzh(mu, mbh, n)
    al = alpha(mu, mbh)
    al2 = al * al
    eps = -8.0 * al2 / (2 * l + 1)
    lp = l + eps
    if lp < 0
        return (0.0, 0.0)
    end
    rp = r_plus(mbh, astar)
    rG = rg(mbh)
    x = sqrt(1.0 - astar * astar)
    y = mu * mu - om0 * om0
    p = -0.5 * (m * astar - 2.0 * rp * om0) / x
    q = 4.0 * om0 * p * rG - 2.0 * (3.0 - x) * al2
    gam_terms = gam_pq_bxzh(p, q, eps, n, l)
    kappab_term = (rG * rG * x * x * y)^(lp + 0.5)
    delta1 = 0.5 * (q / eps - eps - p * 2.0im) * kappab_term * gam_terms
    om = om0 + (eps + delta1) * om1
    return (real(om), imag(om))
end

function h_seidel(l::Integer, m::Integer=1, s::Integer=0)
    num = l * (l * l - m * m)
    denom = 2.0 * (l - 0.5) * (l + 0.5)
    if s > 0
        mabs = abs(m)
        s1 = max(mabs, s)
        s2 = m * s / max(mabs, s)
        l2 = l * l
        num = (l2 - s1 * s1) * (l2 - s * s) * (l2 - s2 * s2)
        denom *= l2 * l
    end
    return num / denom
end

flm_seidel_2(l::Integer, m::Integer=1, s::Integer=0) = h_seidel(l + 1, m, s) - h_seidel(l, m, s) - 1.0

function flm_seidel_4(l::Integer, m::Integer=1, s::Integer=0)
    hl = h_seidel(l, m, s)
    hlp1 = h_seidel(l + 1, m, s)
    hlp2 = h_seidel(l + 2, m, s)
    twol = 2 * l
    l2 = l * l
    lm1 = l - 1
    lp1 = l + 1
    lp2 = l + 2
    res = (hlp1 - lp2 * hlp2 / (twol + 3)) * hlp1 / (2 * lp1)
    res += (hlp1 / lp1 - hl) * hl / twol
    res += lm1 * h_seidel(l - 1, m, s) * hl / (twol * (twol - 1))
    if s > 0
        lp2sq = lp2 * lp2
        lm1sq = lm1 * lm1
        lp1sq = lp1 * lp1
        res += 4.0 * (hlp1 / (lp1sq * lp2sq) - hl / (l2 * lm1sq)) * m * m * s^4 / (l2 * lp1sq)
    end
    return res
end

function alm_approx(c::Number, l::Integer, m::Integer=1, s::Integer=0)
    fvals = [l * (l + 1), flm_seidel_2(l, m, s), flm_seidel_4(l, m, s)]
    expansion = [f * c^(2 * (i - 1)) for (i, f) in enumerate(fvals)]
    return sum(expansion)
end

_calF(s, l, m) = ((0 == s) && (0 == l + 1)) ? 0.0 :
    sqrt(((l + 1)^2 - m * m) / (2 * l + 3) / (2 * l + 1)) *
    sqrt(((l + 1)^2 - s * s) / (l + 1)^2)

_calG(s, l, m) = (0 == l) ? 0.0 :
    sqrt((l * l - m * m) / (4 * l * l - 1)) * sqrt(1.0 - s * s / l / l)

_calH(s, l, m) = ((0 == l) || (0 == s)) ? 0.0 : -m * s / l / (l + 1)
_calA(s, l, m) = _calF(s, l, m) * _calF(s, l + 1, m)
_calD(s, l, m) = _calF(s, l, m) * (_calH(s, l + 1, m) + _calH(s, l, m))
_calB(s, l, m) = _calF(s, l, m) * _calG(s, l + 1, m) + _calG(s, l, m) * _calF(s, l - 1, m) + _calH(s, l, m)^2
_calE(s, l, m) = _calG(s, l, m) * (_calH(s, l - 1, m) + _calH(s, l, m))
_calC(s, l, m) = _calG(s, l, m) * _calG(s, l - 1, m)

swsphericalh_A(s, l, m) = l * (l + 1) - s * (s + 1)

function M_matrix_elem(s, c, m, l, lprime)
    if lprime == l - 2
        return -c * c * _calA(s, lprime, m)
    end
    if lprime == l - 1
        return -c * c * _calD(s, lprime, m) + 2.0 * c * s * _calF(s, lprime, m)
    end
    if lprime == l
        return swsphericalh_A(s, lprime, m) - c * c * _calB(s, lprime, m) + 2.0 * c * s * _calH(s, lprime, m)
    end
    if lprime == l + 1
        return -c * c * _calE(s, lprime, m) + 2.0 * c * s * _calG(s, lprime, m)
    end
    if lprime == l + 2
        return -c * c * _calC(s, lprime, m)
    end
    return zero(complex(c))
end

l_min(s, m) = max(abs(s), abs(m))
ells(s, m, l_max) = collect(l_min(s, m):l_max)

function M_matrix(s, c, m, l_max)
    _ells = ells(s, m, l_max)
    M = Matrix{ComplexF64}(undef, length(_ells), length(_ells))
    for i in eachindex(_ells)
        for j in eachindex(_ells)
            M[i, j] = M_matrix_elem(s, c, m, _ells[i], _ells[j])
        end
    end
    return M
end

sep_consts(s, c, m, l_max) = eigvals(M_matrix(s, c, m, l_max))

function angular_ev(omega::Number, mbh::Real, astar::Real, mu::Real, l::Integer, m::Integer)
    c = rg(mbh) * astar * sqrt(complex(omega * omega - mu * mu))
    if abs(c) > 3
        vals = sort(sep_consts(0, c, m, l); by=real)
        return vals[end]
    end
    return alm_approx(c, l, m, 0)
end

function cfunctions(omega::Number, mbh::Real, astar::Real, mu::Real, alm::Number, m::Integer)
    a = astar
    a2 = a * a
    b = sqrt(1.0 - a2)
    om = rg(mbh) * omega
    om2 = om * om
    mu_r = rg(mbh) * mu
    mu2 = mu_r * mu_r
    q = sqrt(complex(mu2 - om2))
    q = -sign(real(q)) * q
    cN1 = 4.0 * b * q
    cN2 = 0.75 + (2.0 * (b + 1.0) * om2 - (2.0 * b + 1.0) * mu2) / q
    q2 = q * q
    x = (om - 0.5 * m * a) / b
    y = 2.0im * (om + x)
    z = om - 1.0im * q
    z2 = z * z / q
    c0 = 1.0 - y
    c1 = -4.0 + 2.0 * y + 4.0 * (b + 1.0) * q - 2.0 * (q2 + om2) / q
    c2 = 3.0 - y - 2.0 * (q2 - om2) / q
    c3 = 2.0im * z * z2 + a2 * q2 + 2.0im * m * a * q + (z2 + 1.0) * (2.0im * x + 2.0 * b * q - 1.0) - alm
    c4 = z2 * z2 + 2.0im * z2 * (om - x)
    return c0, c1, c2, c3, c4, cN1, cN2
end

alpha_n(n::Integer, c0) = n * n + (c0 + 1.0) * n + c0
beta_n(n::Integer, d1, d3) = -2.0 * n * n + (d1 + 2.0) * n + d3
gamma_n(n::Integer, c2, c4) = n * n + (c2 - 3.0) * n + c4

function continued_fraction(omega::Number, mbh::Real, astar::Real, mu::Real, alm::Number, m::Integer, nmax::Integer=2000)
    c0, c1, c2, c3, c4, cN1, cN2 = cfunctions(omega, mbh, astar, mu, alm, m)
    fr = (-1.0 + 0.0im) + cN1 / sqrt(nmax) + cN2 / nmax
    fr0 = beta_n(0, c1, c3) / alpha_n(0, c0)
    for i in reverse(1:nmax)
        alph = alpha_n(i, c0)
        beta = beta_n(i, c1, c3)
        gam = gamma_n(i, c2, c4)
        fr = gam / (beta - alph * fr)
    end
    return fr0 / fr - 1.0
end

function root_equation(omega::Number, mbh::Real, astar::Real, mu::Real, l::Integer, m::Integer)
    alm = angular_ev(omega, mbh, astar, mu, l, m)
    z = continued_fraction(omega, mbh, astar, mu, alm, m)
    return log(z + 1.0)
end

function _golden_minimize(f, a, b; tol=1e-12, maxiter=160)
    gr = (sqrt(5.0) - 1.0) / 2.0
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = f(c)
    fd = f(d)
    for _ in 1:maxiter
        if abs(b - a) < tol * max(1.0, abs(c), abs(d))
            break
        end
        if fc < fd
            b = d
            d = c
            fd = fc
            c = b - gr * (b - a)
            fc = f(c)
        else
            a = c
            c = d
            fc = fd
            d = a + gr * (b - a)
            fd = f(d)
        end
    end
    x = 0.5 * (a + b)
    return x, f(x)
end

function _coordinate_minimize_2d(cost, x0, y0, xb, yb; passes=8)
    x = clamp(x0, xb[1], xb[2])
    y = clamp(y0, yb[1], yb[2])
    for _ in 1:passes
        x, = _golden_minimize(xx -> cost(xx, y), xb[1], xb[2]; tol=1e-11, maxiter=90)
        y, = _golden_minimize(yy -> cost(x, yy), yb[1], yb[2]; tol=1e-11, maxiter=90)
    end
    return x, y
end

function find_cf_root(mbh::Real, astar::Real, mu::Real, n::Integer=2, l::Integer=1, m::Integer=1; verbose::Bool=false)
    alph = alpha(mu, mbh)
    omR = omegaHyperfine(mu, mbh, astar, n, l, m)
    _, omI = omega_nlm_bxzh(mu, mbh, astar, n, l, m)
    if omR > 0 && omI > 0 && alph > 0
        factor = 0.5 * alph * alph
        om0 = mu * (1.0 - 0.5 * factor * (1.0 / ((n - 1) * (n - 1)) + 1.0 / (n * n)))
        om1 = mu * (1.0 - factor / (n * n))
        if om0 > om1
            om0, om1 = om1, om0
        end

        cost_oR(x) = log(abs(root_equation(x + 1.0im * omI, mbh, astar, mu, l, m)))
        x_best, = _golden_minimize(cost_oR, om0, om1; tol=1e-11, maxiter=120)

        cost(x, lgy) = abs(root_equation(x + 1.0im * 10.0^lgy, mbh, astar, mu, l, m))
        y0 = log10(omI)
        yb = (log10(0.7 * omI), min(log10(10.0 * omI), log10(0.1 * omR)))
        x, lgy = _coordinate_minimize_2d(cost, x_best, y0, (om0, om1), yb; passes=6)
        return x + 1.0im * 10.0^lgy
    else
        if verbose
            println("Estimates for real or imaginary part of omega are not positive")
        end
        return omR + 1.0im * omI
    end
end

function cfm_bhsr_rates(mbh0, mu_max)
    muvals = collect(range(0.0, mu_max, length=250))
    alphvals = alpha.(muvals, mbh0)
    astar0_vals = (4.0 .* alphvals) ./ (1.0 .+ (4.0 .* alphvals.^2))
    roots = [find_cf_root(mbh0, astar, mu, 2, 1, -1) for (mu, astar) in zip(muvals, astar0_vals)]
    return abs.([imag(z) for z in roots] .* rg(mbh0))
end

# =======================================
# Binary Transition Waveforms
# =======================================

function a_tilde_crit(m_i, alpha_in)
    m2 = m_i^2
    return (4.0 * m_i * alpha_in) / (m2 + 4.0 * alpha_in^2)
end

function q_c(alpha_in, m_i)
    m = float(m_i)
    one_minus = 1.0 - alpha_in / m
    inner = 1.0 - (4.0 * alpha_in / m * one_minus)^2
    denom = m^2 * (1.0 - sqrt(inner))
    if denom == 0.0
        return 0.0
    end
    return 8.0 * alpha_in^2 * one_minus / denom - 1.0
end

function h0_from_params(qc, M, r, alpha_in, Omega0)
    return 24.0 * G * (qc * M / r) * alpha_in^-4 * (G * M * Omega0)^2
end

function gamma_rate(q, M, Omega0)
    return Omega0^2 * (96.0 / 5.0) * (q / (q + 1.0)^(1.0 / 3.0)) * (G * M * Omega0)^(5.0 / 3.0)
end

z_parameter(eta, Delta_m, gamma) = eta^2 / (abs(Delta_m) * gamma)

function z_scaling_211_to_21m1(alpha_in, q)
    return 7.0 * (1.81 / (1.0 + 4.0 * alpha_in^2))^(1.0 / 3.0) *
           q^(1.0 / 3.0) *
           (2.0 / (1.0 + q))^(5.0 / 3.0) *
           (0.45 / alpha_in)^(11.0 / 3.0)
end

fc_from_Omega0(Omega0) = (2.0 / pi) * Omega0

function psi_plus(f, r, f0, Delta_m, gamma)
    return f .* r .+ ((f .- f0).^2) ./ (4.0 * abs(Delta_m) * gamma) .- pi / 4.0
end

function htilde_plus(
    f,
    M,
    r,
    alpha_in,
    Omega0,
    q,
    m_i,
    m_f,
    eta,
    Gamma_abs;
    use_z_scaling=false,
    numerical_qc=true,
)
    f = Float64.(collect(f))
    Delta_m = abs(m_f - m_i)

    acrit = a_tilde_crit(m_i, alpha_in)
    gamma = gamma_rate(q, M, Omega0)
    boson_mass = alpha_in / (G * M)

    z = use_z_scaling ? z_scaling_211_to_21m1(alpha_in, q) : z_parameter(eta, Delta_m, gamma)

    qc = if numerical_qc
        compute_cloud_mass_numerical(alpha_in, M, boson_mass; spin=0.99) / M
    else
        q_c(alpha_in, m_i)
    end

    h0 = h0_from_params(qc, M, r, alpha_in, Omega0)

    f_c = fc_from_Omega0(Omega0)
    f0 = f_c
    phase = psi_plus(f, r, f0, Delta_m, gamma)

    denom = sqrt(z) ./ (abs(Gamma_abs) .- 1.0im .* pi .* (f .- f_c))
    envelope = exp(-pi * z) .* exp.(-2.0 * z .* atan.(pi .* (f .- f_c) ./ abs(Gamma_abs)))

    function with_inclination(iota)
        pref = h0 * (1.0 + cos(iota)^2) * sqrt(pi) * Delta_m^2
        return abs.(pref .* 1.0im .* exp.(1.0im .* phase) .* envelope .* denom) ./ 2.417987242e14
    end

    return with_inclination
end

end # module

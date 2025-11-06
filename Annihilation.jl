


using SpecialFunctions, Test
import SpecialFunctions: expinti


"""
    expinti(z::Complex)

Exponential integral expinti(z) for complex argument z, matching Python's mpmath.ei(z)
"""
function expinti(z::Complex)
    E1 = expint(1, -z)
    if imag(z) > 0
        return -E1 + im * π
    elseif imag(z) < 0
        return -E1 - im * π
    else
        # z real: approach from above
        return -E1
    end
end


using SpecialFunctions, CairoMakie

# Requires omega in Hz, M in solar masses, alpha = dimensionless, r in kpc
# r is distance from us to binary

# H_ann gives the change in frequency from \omega_ann. Centered in the MHz, but with \delta f being extremely small.

const Mp = 1.220890e19
const G = 1 / Mp^2


# Define called on functions
function omega_ann(mua, alpha, n)
    return 2 * mua * (1 - alpha^2 / (2 * n^2))
end

function gamma_a(l, alpha, M)
    if l == 1
        p = 17
    else
        p = 4 * l + 1
    end
    rg = G * M
    return G * 1e-10 / rg^3 * (((alpha / l) * 0.5)^p + ((alpha / l) * 0.5)^(p + 1))
end

function h_ann(omega, M, n, l, alpha, r)

    # Define Conversions and Constants
    Mp_local = 1.220890e19
    G_local = 1 / Mp_local^2
    M_sun = 2e30 * 5.61e26 # Solar mass -> GeV
    Mpc_to_Gev = 1.56e38        # 1 Mpc in GeV
    kpc_to_Gev = Mpc_to_Gev / 1000
    Hz_to_Gev = 4.1357e-24 # Hz -> GeV
    hbar_GeV = 6.582e-25
    
    # Convert to natural units
    M_ann = M * M_sun # Solar mass -> GeV
    mua = alpha / (G_local * M_ann) # Mu in GeV
    r_kpc = r * kpc_to_Gev # Kpc -> GeV

    omega_a = omega_ann(mua, alpha, n)
    Gamma = gamma_a(l, alpha, M_ann)

    # Approximate maximum population of states
    N_max = 10.0^76 * (M/10.0)^2
    omega_gev = omega .* Hz_to_Gev
    z = 1im * omega_gev / (Gamma * N_max)  # argument for Ei and exponential
    pref = 1/(2*π) * sqrt(4.0 * G_local / (Gamma * r_kpc^2 * omega_a))
    Ei = -2.0 * expinti.(z) - log.(-1im * omega_gev) + log.(1im * omega_gev) + 1im * π * sign.(real.(omega_gev))

    freq = omega_gev / (2*π*Hz_to_Gev)
    
    # frequency strain has units of h*GeV, need to convert to Hz and then divide by timescale of event (2e6) to match discrete fft
    return pref * exp.(z) .* Ei * hbar_GeV, freq
end




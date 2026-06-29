include(joinpath(@__DIR__, "..", "functions.jl"))

using DelimitedFiles
using .HenryWF

const TEST_DIR = joinpath(@__DIR__, "test_files", "julia")

function save_real(path, values)
    array = values isa Number ? [Float64(values)] : Float64.(collect(values))
    writedlm(path, reshape(array, :, 1))
end

function save_complex(path, values)
    array = values isa Number ? [ComplexF64(values)] : ComplexF64.(collect(values))
    writedlm(path, hcat(real.(array), imag.(array)))
end

function main()
    mkpath(TEST_DIR)

    alpha_in = 0.42
    q = 0.4
    M_sol = 10.0
    M = M_sol * Msol_in_eV
    r = 1.0e28
    m_i = 1
    m_f = -1
    n = 2
    l_i = 1
    Omega0 = Omega0_binary_natural_unit(m_i, alpha_in, G, M, n, l_i)
    eta = eta_parameter(alpha_in, q, M_sol)
    Gamma_abs = 1.0e-20
    f = collect(range(0.5, 1.5, length=32)) .* fc_from_Omega0(Omega0)
    h_of_iota = htilde_plus(
        f, M, r, alpha_in, Omega0, q, m_i, m_f, eta, Gamma_abs;
        use_z_scaling=false, numerical_qc=false,
    )

    save_real(joinpath(TEST_DIR, "eta.txt"), eta)
    save_real(joinpath(TEST_DIR, "omega0_binary.txt"), Omega0)
    save_real(joinpath(TEST_DIR, "gamma_rate.txt"), gamma_rate(q, M, Omega0))
    save_real(joinpath(TEST_DIR, "q_c.txt"), q_c(alpha_in, m_i))
    save_real(joinpath(TEST_DIR, "z_scaling.txt"), z_scaling_211_to_21m1(alpha_in, q))
    save_real(joinpath(TEST_DIR, "f_grid.txt"), f)
    save_real(joinpath(TEST_DIR, "htilde_iota0.txt"), h_of_iota(0.0))
    save_real(joinpath(TEST_DIR, "htilde_iota1.txt"), h_of_iota(1.0))

    mu = 1.0e-12
    mbh = 1.0e-6
    astar = 0.7
    omega = 0.99 * mu + 1.0e-20im
    alm = angular_ev(omega, mbh, astar, mu, 1, 1)
    save_real(joinpath(TEST_DIR, "rg.txt"), rg(mbh))
    save_real(joinpath(TEST_DIR, "alpha_bhsr.txt"), alpha(mu, mbh))
    save_real(joinpath(TEST_DIR, "r_plus.txt"), r_plus(mbh, astar))
    save_real(joinpath(TEST_DIR, "omega_hyperfine.txt"), omegaHyperfine(mu, mbh, astar, 2, 1, 1))
    save_real(joinpath(TEST_DIR, "omega0_bxzh.txt"), omega0_bxzh(mu, mbh, 2))
    save_real(joinpath(TEST_DIR, "omega1_bxzh.txt"), omega1_bxzh(mu, mbh, 2))
    save_complex(joinpath(TEST_DIR, "omega_nlm_bxzh.txt"), [complex(omega_nlm_bxzh(mu, mbh, astar, 2, 1, 1)...)])
    save_complex(joinpath(TEST_DIR, "angular_ev.txt"), [alm])
    save_complex(joinpath(TEST_DIR, "continued_fraction.txt"), [continued_fraction(omega, mbh, astar, mu, alm, 1, 80)])
end

main()

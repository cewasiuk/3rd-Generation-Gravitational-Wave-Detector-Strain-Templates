include(joinpath(@__DIR__, "..", "annihilation.jl"))
include(joinpath(@__DIR__, "..", "level_transition.jl"))

using DelimitedFiles
using .ChrisWF

const TEST_DIR = joinpath(@__DIR__, "test_files", "julia")

function save_real(path, values)
    array = values isa Number ? [Float64(values)] : Float64.(collect(values))
    writedlm(path, reshape(array, :, 1))
end

function save_complex(path, values)
    array = values isa Number ? [ComplexF64(values)] : ComplexF64.(collect(values))
    writedlm(path, hcat(real.(array), imag.(array)))
end

function save_row(path, values)
    array = values isa Number ? [Float64(values)] : Float64.(collect(values))
    writedlm(path, reshape(array, 1, :))
end

function main()
    mkpath(TEST_DIR)

    level = iso_gatom_level_tr_strain(
        M_solar=1e-6,
        a_spin=0.999999,
        alpha=0.75,
        ne=6,
        ng=5,
        distance_kpc=1.0,
        N_g0=1e-6,
        n_time=512,
        n_fft=4096,
        n_top=120,
        verbose=false,
    )

    annihilation = iso_gatom_ann_strain(
        M_solar=3.1e-4,
        mua=2e-16,
        n=4,
        l=nothing,
        alpha=nothing,
        distance_kpc=1.0,
        iota=0.0,
        phase=0.0,
        f_min_Hz=1e7,
        f_max_Hz=1e9,
        n_f=1024,
        verbose=false,
    )

    for (name, values) in level
        path = joinpath(TEST_DIR, "level_$(name).txt")
        if name == "lorentz_params"
            save_row(path, values)
        elseif values isa AbstractArray{<:Complex}
            save_complex(path, values)
        else
            save_real(path, values)
        end
    end

    for (name, values) in annihilation
        path = joinpath(TEST_DIR, "ann_$(name).txt")
        if values isa Number
            save_real(path, [values])
        elseif values isa AbstractArray{<:Complex}
            save_complex(path, values)
        else
            save_real(path, values)
        end
    end
end

main()

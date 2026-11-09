"""
    FQHTranslator <: AbstractTranslator

Translator for Fractional Quantum Hall states on a cylinder.
Translates quantum numbers (sectors) across unit cell boundaries.

In the Landau gauge on a cylinder, translating by `shift` unit cells
shifts the momentum quantum number by:
    ΔP = shift * L_cell * N_particles
where `L_cell` (norbitals) is the number of orbitals in the unit cell, and `N_particles` is the
particle number charge.
"""
struct FQHTranslator <: AbstractTranslator
    norbitals::Int
    nparticles::Int
    trs::Bool
end

FQHTranslator(norbitals::Int; nparticles::Int = 0, trs::Bool = false) =
    FQHTranslator(norbitals, nparticles, trs)

function translate_sector(c::ProductSector{T}, tr::FQHTranslator, shift::Integer, dual::Bool = false) where {T <: Tuple}
    # c[1] is charge N (U1Irrep)
    # c[2] is momentum P (U1Irrep)
    N_ch = c[1].charge
    P_ch = c[2].charge
    sign_val = tr.trs ? -1 : 1
    new_N = N_ch + shift * tr.nparticles
    new_P = P_ch + shift * tr.norbitals * N_ch * sign_val
    if length(c.sectors) == 2
        return U1Irrep(new_N) ⊠ U1Irrep(new_P)
    else
        other_sectors = c.sectors[3:end]
        return ProductSector((U1Irrep(new_N), U1Irrep(new_P), other_sectors...))
    end
end

function translate_sector(c::U1Irrep, tr::FQHTranslator, shift::Integer, dual::Bool = false)
    return U1Irrep(c.charge + shift * tr.norbitals)
end

function translate_sector(c::TensorKit.Trivial, tr::FQHTranslator, shift::Integer, dual::Bool = false)
    return c
end

translate_sector(c, tr::FQHTranslator, shift::Integer) = translate_sector(c, tr, shift, false)

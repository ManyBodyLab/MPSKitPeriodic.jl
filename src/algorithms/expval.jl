function MPSKit.expectation_value(
        ψ::InfinitePeriodicMPS, H::InfinitePeriodicMPOHamiltonian,
        envs::AbstractMPSEnvironments = environments(ψ, H)
    )
    return sum(1:length(ψ)) do site
        return MPSKit.contract_mpo_expval(
            ψ.AC[site], envs.GLs[site], H[site][:, 1, 1, end], envs.GRs[site][end]
        )
    end
end

# MPO tensor
#-----------
function MPSKit.expectation_value(
        ψ::InfinitePeriodicMPS, mpo_pair::Tuple{InfinitePeriodicMPO, Pair{Int, <:MPOTensor}},
        envs::AbstractMPSEnvironments = environments(ψ, first(mpo_pair))
    )
    O_mpo, (site, O) = mpo_pair
    return MPSKit.contract_mpo_expval(ψ.AC[site], envs.GLs[site], O, envs.GRs[site]) /
        MPSKit.expectation_value(ψ, O_mpo, envs)
end

# DenseMPO
# --------
function MPSKit.expectation_value(ψ::InfinitePeriodicMPS, mpo::InfinitePeriodicMPO, envs...)
    return MPSKit.expectation_value(convert(MultilineMPS, ψ), convert(MultilineMPO, mpo), envs...)
end

# Density matrices
# ----------------
function MPSKit.expectation_value(ρ::InfinitePeriodicMPO, args...)
    return MPSKit.expectation_value(convert(InfinitePeriodicMPS, ρ), args...)
end

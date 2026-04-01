MPSKit.entropy(state::InfinitePeriodicMPS) = map(Base.Fix1(MPSKit.entropy, state), 1:length(state))
function MPSKit.entropy(state::InfinitePeriodicMPS, loc::Int)
    return MPSKit.entropy(entanglement_spectrum(state, loc))
end

function MPSKit.entanglement_spectrum(st::InfinitePeriodicMPS, site::Int = 0)
    checkbounds(st, site)
    return LinearAlgebra.svdvals(st.C[site])
end

function MPSKit.marek_gap(above::InfinitePeriodicMPS; tol_angle = 0.1, kwargs...)
    spectrum = MPSKit.transfer_spectrum(above; kwargs...)
    return MPSKit.marek_gap(spectrum; tol_angle)
end

function MPSKit.correlation_length(above::InfinitePeriodicMPS; kwargs...)
    ϵ, = MPSKit.marek_gap(above; kwargs...)
    return 1 / ϵ
end

function MPSKit.variance(
        state::InfinitePeriodicMPS, H::InfinitePeriodicMPOHamiltonian, envs = environments(state, H)
    )
    e_local = map(1:length(state)) do i
        return contract_mpo_expval(
            state.AC[i], envs.GLs[i], H[i][:, :, :, end], envs.GRs[i][end]
        )
    end
    lattice = physicalspace(state)
    ## TODO: Make sure this function calls the correct translation function
    H_renormalized = InfinitePeriodicMPOHamiltonian(
        lattice, i => e * id(storagetype(eltype(H)), lattice[i]) for (i, e) in enumerate(e_local)
    )

    return real(MPSKit.expectation_value(state, (H - H_renormalized)^2))
end

function MPSKit.transfer_spectrum(
        above::InfinitePeriodicMPS; below = above, tol = MPSKit.Defaults.tol, num_vals = 20,
        sector = leftunit(above)
    )
    init = MPSKit.randomize!(
        similar(
            above.AL[1], left_virtualspace(below, 1),
            spacetype(above)(sector => 1)' * left_virtualspace(above, 1)
        )
    )

    transferspace = fuse(left_virtualspace(above, 1) * left_virtualspace(below, 1)')
    num_vals = min(dim(transferspace, sector), num_vals) # we can ask at most this many values
    eigenvals, eigenvecs, convhist = MPSKit.eigsolve(
        flip(MPSKit.TransferMatrix(above.AL, below.AL)), init, num_vals, :LM; tol = tol
    )
    convhist.converged < num_vals &&
        @warn "correlation length failed to converge: normres = $(convhist.normres)"

    return eigenvals
end

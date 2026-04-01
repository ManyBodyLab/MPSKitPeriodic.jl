using MPSKit.Defaults

function MPSKit.environment_alg(
        ::InfinitePeriodicMPS, ::InfinitePeriodicMPO, ::InfinitePeriodicMPS;
        tol = Defaults.tol, maxiter = Defaults.maxiter, krylovdim = Defaults.krylovdim,
        verbosity = Defaults.VERBOSE_NONE, eager = true
    )
    return MPSKit.Arnoldi(; tol, maxiter, krylovdim, verbosity, eager)
end
function MPSKit.environment_alg(
        below, ::InfinitePeriodicMPOHamiltonian, above;
        tol = Defaults.tol, maxiter = Defaults.maxiter, krylovdim = Defaults.krylovdim,
        verbosity = Defaults.VERBOSE_NONE
    )
    max_krylovdim = ceil(Int, dim(left_virtualspace(above, 1)) * dim(left_virtualspace(below, 1)))
    return MPSKit.GMRES(; tol, maxiter, krylovdim = min(max_krylovdim, krylovdim), verbosity)
end

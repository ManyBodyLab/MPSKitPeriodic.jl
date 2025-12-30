using Literate: Literate
using MPSKitPeriodic

Literate.markdown(
    joinpath(pkgdir(MPSKitPeriodic), "docs", "files", "README.jl"),
    joinpath(pkgdir(MPSKitPeriodic));
    flavor = Literate.CommonMarkFlavor(),
    name = "README",
)

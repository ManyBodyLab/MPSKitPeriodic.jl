using Literate: Literate
using MPSKitPeriodic

Literate.markdown(
    joinpath(pkgdir(MPSKitPeriodic), "docs", "files", "README.jl"),
    joinpath(pkgdir(MPSKitPeriodic), "docs", "src");
    flavor = Literate.DocumenterFlavor(),
    name = "index",
)

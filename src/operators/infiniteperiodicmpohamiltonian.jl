const InfinitePeriodicMPOHamiltonian{O <: MPOTensor, F} = MPOHamiltonian{O, PeriodicVector{O, F}}

MPSKit.GeometryStyle(::Type{<:InfinitePeriodicMPOHamiltonian}) = InfiniteChainStyle()

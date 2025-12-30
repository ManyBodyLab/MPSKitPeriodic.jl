const InfinitePeriodicMPO{O, F} = MPO{O, PeriodicVector{O, F}}

MPSKit.GeometryStyle(::Type{<:InfinitePeriodicMPO}) = InfiniteChainStyle()

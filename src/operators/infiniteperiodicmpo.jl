const InfinitePeriodicMPO{O, F, G} = MPO{O, PeriodicVector{O, F, G}}

Base.isfinite(::Type{<:InfinitePeriodicMPO}) = false
MPSKit.GeometryStyle(::Type{<:InfinitePeriodicMPO}) = InfiniteChainStyle()

function InfinitePeriodicMPO(Ws::AbstractVector)
    return InfinitePeriodicMPO(PeriodicVector(Ws))
end
function InfinitePeriodicMPO(Os::PeriodicVector{O, F, G}) where {O, F, G}
    for i in eachindex(Os)
        right_virtualspace(Os[i]) == left_virtualspace(Os[mod1(i + 1, end)]) ||
            throw(SpaceMismatch("unmatching virtual spaces at site $i"))
    end
    return InfinitePeriodicMPO{O, F, G}(Os)
end

function InfinitePeriodicMPO(H::MPSKit.InfiniteMPO)
    return InfinitePeriodicMPO(H.O)
end

const InfinitePeriodicMPO{O, F} = MPO{O, PeriodicVector{O, F}}

Base.isfinite(::Type{<:InfinitePeriodicMPO}) = false
MPSKit.GeometryStyle(::Type{<:InfinitePeriodicMPO}) = InfiniteChainStyle()

function InfinitePeriodicMPO(Ws::AbstractVector)
    return InfinitePeriodicMPO(PeriodicVector(Ws))
end
function InfinitePeriodicMPO(Os::PeriodicVector{O, F}) where {O, F}
    for i in eachindex(Os)
        right_virtualspace(Os[i]) == left_virtualspace(Os[mod1(i + 1, end)]) ||
            throw(SpaceMismatch("unmatching virtual spaces at site $i"))
    end
    return InfinitePeriodicMPO{O, F}(Os)
end

function InfinitePeriodicMPO(H::MPSKit.InfiniteMPO)
    return InfinitePeriodicMPO(H.O)
end


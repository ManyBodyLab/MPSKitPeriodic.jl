const InfinitePeriodicMPOHamiltonian{O <: MPOTensor, F, G} = MPOHamiltonian{O, PeriodicVector{O, F, G}}

Base.isfinite(::Type{<:InfinitePeriodicMPOHamiltonian}) = false
MPSKit.GeometryStyle(::Type{<:InfinitePeriodicMPOHamiltonian}) = InfiniteChainStyle()

function InfinitePeriodicMPOHamiltonian(Ws::AbstractVector, translator=x->x[1])
    return InfinitePeriodicMPOHamiltonian(PeriodicVector(Ws, translator))
end
function InfinitePeriodicMPOHamiltonian(Ws::PeriodicVector{O, F, G}, translator=nothing) where {O <: MPOTensor, F, G}
    for i in eachindex(Ws)
        right_virtualspace(Ws[i]) == left_virtualspace(Ws[i + 1]) ||
            throw(ArgumentError("The virtual spaces of the MPO tensors at site $i do not match."))
    end
    return InfinitePeriodicMPOHamiltonian{O, F, G}(Ws)
end


function MPSKit.isidentitylevel(H::InfinitePeriodicMPOHamiltonian{<:JordanMPOTensor}, i::Int)
    if i == 1 || i == size(H[1], 1)
        return true
    else
        return all(H.A) do A
            return haskey(A, CartesianIndex(i - 1, 1, 1, i - 1)) &&
                A[i - 1, 1, 1, i - 1] isa BraidingTensor
        end
    end
end
function MPSKit.isemptylevel(H::InfinitePeriodicMPOHamiltonian, i::Int)
    return any(parent(H)) do h
        return !haskey(h, CartesianIndex(i, 1, 1, i))
    end
end

function InfinitePeriodicMPOHamiltonian(H::MPSKit.InfiniteMPOHamiltonian)
    return InfinitePeriodicMPOHamiltonian(H.W)
end

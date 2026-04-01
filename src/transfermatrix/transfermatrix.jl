#the product of transfer matrices is its own type
struct PeriodicProductTransferMatrix{T <: AbstractTransferMatrix, P <: PeriodicVector{T}} <: AbstractTransferMatrix
    tms::P
end

PeriodicProductTransferMatrix(v::AbstractVector) = PeriodicProductTransferMatrix(convert(PeriodicVector, v));

# flip em
TensorKit.flip(tm::PeriodicProductTransferMatrix) = PeriodicProductTransferMatrix(flip.(reverse(tm.tms)));

# TransferMatrix acting as a function
function (d::PeriodicProductTransferMatrix)(vec) 
    return d.tms.map(foldr((a, b) -> a(b), d.tms; init = vec), 1)
end

function MPSKit.TransferMatrix(a::PeriodicVector, b, c::PeriodicVector, isflipped = false)
    tot = PeriodicProductTransferMatrix(PeriodicVector(MPSKit.TransferMatrix.(a, b, c), a.map))
    return isflipped ? flip(tot) : tot
end


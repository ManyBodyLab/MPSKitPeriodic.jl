function MPSKit.makefullrank!(A::PeriodicVector{<:GenericMPSTensor}; alg_leftorth = MPSKit.Defaults.alg_qr(), alg_rightorth = MPSKit.Defaults.alg_lq())
    while true
        i = findfirst(!isfullrank, A)
        isnothing(i) && break
        if !isfullrank(A[i]; side = :left)
            L, Q = right_orth!(_transpose_tail(A[i]); alg = alg_rightorth)
            A[i] = _transpose_front(Q)
            A[i - 1] = A[i - 1] * L
        else
            A[i], R = left_orth!(A[i]; alg = alg_leftorth)
            A[i + 1] = _transpose_front(R * _transpose_tail(A[i + 1]))
        end
    end
    return A
end

function MPSKit.makefullrank!(virtualspaces::PeriodicVector{S}, physicalspaces::PeriodicVector{S}) where {S <: ElementarySpace}
    haschanged = true
    while haschanged
        haschanged = false
        for i in 1:length(virtualspaces)
            Vmax = fuse(virtualspaces[i - 1], physicalspaces[i - 1])
            if !(virtualspaces[i] ≾ Vmax)
                virtualspaces[i] = infimum(virtualspaces[i], Vmax)
                haschanged = true
            end
        end
        for i in reverse(1:length(virtualspaces))
            Vmax = fuse(dual(physicalspaces[i - 1]), virtualspaces[i])
            if !(virtualspaces[i - 1] ≾ Vmax)
                virtualspaces[i - 1] = infimum(virtualspaces[i - 1], Vmax)
                haschanged = true
            end
        end
    end

    return virtualspaces
end

function MPSKit.gaugefix!(
        ψ::InfinitePeriodicMPS, A, C₀ = ψ.C[end];
        order = :LR, timeroutput = MPSKit.NoTimerOutput(),
        backend::MPSKit.AbstractBackend = MPSKit.DefaultBackend(), allocator = MPSKit.DefaultAllocator(),
        kwargs...
    )
    alg = if order === :LR || order === :RL
        MPSKit.MixedCanonical(; order, kwargs...)
    elseif order === :L
        MPSKit.LeftCanonical(; kwargs...)
    elseif order === :R
        MPSKit.RightCanonical(; kwargs...)
    else
        throw(ArgumentError("Invalid order: $order"))
    end

    return MPSKit.gaugefix!(ψ, A, C₀, alg; timeroutput, backend, allocator)
end

# expert mode: actual implementation
function MPSKit.gaugefix!(
        ψ::InfinitePeriodicMPS, A, C₀, alg::MPSKit.MixedCanonical;
        timeroutput = MPSKit.NoTimerOutput(),
        backend::MPSKit.AbstractBackend = MPSKit.DefaultBackend(), allocator = MPSKit.DefaultAllocator()
    )
    if alg.order === :LR
        MPSKit.gaugefix!(ψ, A, C₀, alg.alg_leftcanonical; timeroutput, backend, allocator)
        MPSKit.gaugefix!(ψ, ψ.AL, ψ.C[end], alg.alg_rightcanonical; timeroutput, backend, allocator)
    elseif alg.order === :RL
        MPSKit.gaugefix!(ψ, A, C₀, alg.alg_rightcanonical; timeroutput, backend, allocator)
        MPSKit.gaugefix!(ψ, ψ.AR, ψ.C[end], alg.alg_leftcanonical; timeroutput, backend, allocator)
    else
        throw(ArgumentError("Invalid order: $(alg.order)"))
    end
    return ψ
end
function MPSKit.gaugefix!(
        ψ::InfinitePeriodicMPS, A, C₀, alg::MPSKit.LeftCanonical;
        timeroutput = MPSKit.NoTimerOutput(),
        backend::MPSKit.AbstractBackend = MPSKit.DefaultBackend(), allocator = MPSKit.DefaultAllocator()
    )
    MPSKit.uniform_leftorth!((ψ.AL, ψ.C), A, C₀, alg; timeroutput, backend, allocator)
    return ψ
end
function MPSKit.gaugefix!(
        ψ::InfinitePeriodicMPS, A, C₀, alg::MPSKit.RightCanonical;
        timeroutput = MPSKit.NoTimerOutput(),
        backend::MPSKit.AbstractBackend = MPSKit.DefaultBackend(), allocator = MPSKit.DefaultAllocator()
    )
    MPSKit.uniform_rightorth!((ψ.AR, ψ.C), A, C₀, alg; timeroutput, backend, allocator)
    return ψ
end

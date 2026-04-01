function MPSKit.gaugefix!(ψ::InfinitePeriodicMPS, A, C₀ = ψ.C[end]; order = :LR, kwargs...)
    alg = if order === :LR || order === :RL
        MPSKit.MixedCanonical(; order, kwargs...)
    elseif order === :L
        MPSKit.LeftCanonical(; kwargs...)
    elseif order === :R
        MPSKit.RightCanonical(; kwargs...)
    else
        throw(ArgumentError("Invalid order: $order"))
    end

    return MPSKit.gaugefix!(ψ, A, C₀, alg)
end

# expert mode: actual implementation
function MPSKit.gaugefix!(ψ::InfinitePeriodicMPS, A, C₀, alg::MPSKit.MixedCanonical)
    if alg.order === :LR
        MPSKit.gaugefix!(ψ, A, C₀, alg.alg_leftcanonical)
        MPSKit.gaugefix!(ψ, ψ.AL, ψ.C[end], alg.alg_rightcanonical)
    elseif alg.order === :RL
        MPSKit.gaugefix!(ψ, A, C₀, alg.alg_rightcanonical)
        MPSKit.gaugefix!(ψ, ψ.AR, ψ.C[end], alg.alg_leftcanonical)
    else
        throw(ArgumentError("Invalid order: $(alg.order)"))
    end
    return ψ
end
function MPSKit.gaugefix!(ψ::InfinitePeriodicMPS, A, C₀, alg::MPSKit.LeftCanonical)
    MPSKit.uniform_leftorth!((ψ.AL, ψ.C), A, C₀, alg)
    return ψ
end
function MPSKit.gaugefix!(ψ::InfinitePeriodicMPS, A, C₀, alg::MPSKit.RightCanonical)
    MPSKit.uniform_rightorth!((ψ.AR, ψ.C), A, C₀, alg)
    return ψ
end

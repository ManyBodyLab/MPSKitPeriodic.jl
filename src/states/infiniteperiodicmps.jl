"""
    InfinitePeriodicMPS{A<:GenericMPSTensor,B<:MPSBondTensor} <: AbtractMPS

Type that represents an infinite Matrix Product State.

## Fields
- `AL` -- left-gauged MPS tensors
- `AR` -- right-gauged MPS tensors
- `AC` -- center-gauged MPS tensors
- `C` -- gauge tensors

## Notes
By convention, we have that:
- `AL[i] * C[i]` = `AC[i]` = `C[i-1] * AR[i]`
- `AL[i]' * AL[i] = 1`
- `AR[i] * AR[i]' = 1`

---

## Constructors
    InfinitePeriodicMPS([f, eltype], physicalspaces::Vector{<:Union{S, CompositeSpace{S}},
                virtualspaces::Vector{<:Union{S, CompositeSpace{S}};
                kwargs...) where {S<:ElementarySpace}
    InfinitePeriodicMPS(As::AbstractVector{<:GenericMPSTensor}; kwargs...)
    InfinitePeriodicMPS(ALs::AbstractVector{<:GenericMPSTensor}, C₀::MPSBondTensor;
                kwargs...)

Construct an MPS via a specification of physical and virtual spaces, or from a list of
tensors `As`, or a list of left-gauged tensors `ALs`.

### Arguments
- `As::AbstractVector{<:GenericMPSTensor}`: vector of site tensors
- `ALs::AbstractVector{<:GenericMPSTensor}`: vector of left-gauged site tensors
- `C₀::MPSBondTensor`: initial gauge tensor

- `f::Function=rand`: initializer function for tensor data
- `eltype::Type{<:Number}=ComplexF64`: scalar type of tensors

- `physicalspaces::AbstractVector{<:Union{S, CompositeSpace{S}}`: list of physical spaces
- `virtualspaces::AbstractVector{<:Union{S, CompositeSpace{S}}`: list of virtual spaces

### Keywords
- `tol`: gauge fixing tolerance
- `maxiter`: gauge fixing maximum iterations
"""
struct InfinitePeriodicMPS{A <: GenericMPSTensor, B <: MPSBondTensor, P1 <: PeriodicVector{A}, P2 <: PeriodicVector{A}, P3 <: PeriodicVector{B}, P4 <: PeriodicVector{A}} <: AbstractMPS
    AL::P1
    AR::P2
    C::P3
    AC::P4
    function InfinitePeriodicMPS{A, B}(
            AL::PeriodicVector{A}, AR::PeriodicVector{A},
            C::PeriodicVector{B}, AC::PeriodicVector{A} = AL .* C
        ) where {A <: GenericMPSTensor, B <: MPSBondTensor}
        # verify lengths are compatible
        L = length(AL)
        L == length(AR) == length(C) == length(AC) ||
            throw(ArgumentError("incompatible lengths of AL, AR, C, and AC"))
        # verify tensors are compatible
        spacetype(A) == spacetype(B) ||
            throw(SpaceMismatch("incompatible space types of AL and C"))
        P1 = typeof(AL)
        P2 = typeof(AR)
        P3 = typeof(C)
        P4 = typeof(AC)
        return new{A, B, P1, P2, P3, P4}(AL, AR, C, AC)
    end
    function InfinitePeriodicMPS(
            AL::PeriodicVector{A}, AR::PeriodicVector{A},
            C::PeriodicVector{B}, AC::PeriodicVector{A} = AL .* C
        ) where {A <: GenericMPSTensor, B <: MPSBondTensor}
        # verify lengths are compatible
        L = length(AL)
        L == length(AR) == length(C) == length(AC) ||
            throw(ArgumentError("incompatible lengths of AL, AR, C, and AC"))
        # verify tensors are compatible
        spacetype(A) == spacetype(B) ||
            throw(SpaceMismatch("incompatible space types of AL and C"))

        for i in 1:L
            N = numind(AL[i])
            N == numind(AR[i]) == numind(AC[i]) ||
                throw(SpaceMismatch("incompatible spaces at site $i"))

            # verify that the physical spaces are compatible
            phys_ind = 2:(N - 1)
            all(
                space.(Ref(AL[i]), phys_ind) .== space.(Ref(AR[i]), phys_ind) .==
                    space.(Ref(AC[i]), phys_ind)
            ) ||
                throw(SpaceMismatch("incompatible physical spaces at site $i"))

            # verify that the virtual spaces are compatible
            space(AL[i], 1) == dual(space(AL[i - 1], N)) &&
                space(AR[i], 1) == dual(space(AR[i - 1], N)) &&
                space(AC[i], 1) == space(AL[i], 1) &&
                space(AC[i], N) == space(AR[i], N) &&
                space(C[i], 1) == dual(space(AL[i], N)) &&
                space(AR[i], 1) == dual(space(C[i - 1], 2)) ||
                throw(SpaceMismatch("incompatible virtual spaces at site $i"))
            # verify that the spaces are non-zero
            dim(space(AL[i])) > 0 && dim(space(C[i])) > 0 ||
                @warn "no fusion channels available at site $i"
        end
        P1 = typeof(AL)
        P2 = typeof(AR)
        P3 = typeof(C)
        P4 = typeof(AC)
        return new{A, B, P1, P2, P3, P4}(AL, AR, C, AC)
    end
end

#===========================================================================================
Constructors
===========================================================================================#

function InfinitePeriodicMPS(
        AL::AbstractVector{A}, AR::AbstractVector{A}, C::AbstractVector{B},
        AC::AbstractVector{A} = AL .* C
    ) where {A <: GenericMPSTensor, B <: MPSBondTensor}
    return InfinitePeriodicMPS(
        PeriodicVector{A}(AL), PeriodicVector{A}(AR),
        PeriodicVector{B}(C), PeriodicVector{A}(AC)
    )
end

function InfinitePeriodicMPS(
        pspaces::AbstractVector{S}, Dspaces::AbstractVector{S};
        kwargs...
    ) where {S <: IndexSpace}
    return InfinitePeriodicMPS(MPSTensor.(pspaces, circshift(Dspaces, 1), Dspaces); kwargs...)
end
function InfinitePeriodicMPS(
        f, elt::Type{<:Number}, pspaces::AbstractVector{S}, Dspaces::AbstractVector{S};
        kwargs...
    ) where {S <: IndexSpace}
    return InfinitePeriodicMPS(
        MPSTensor.(f, elt, pspaces, circshift(Dspaces, 1), Dspaces);
        kwargs...
    )
end
InfinitePeriodicMPS(d::S, D::S) where {S <: Union{Int, <:IndexSpace}} = InfinitePeriodicMPS([d], [D])
function InfinitePeriodicMPS(
        f, elt::Type{<:Number}, d::S, D::S
    ) where {S <: Union{Int, <:IndexSpace}}
    return InfinitePeriodicMPS(f, elt, [d], [D])
end
function InfinitePeriodicMPS(ds::AbstractVector{Int}, Ds::AbstractVector{Int})
    return InfinitePeriodicMPS(ComplexSpace.(ds), ComplexSpace.(Ds))
end
function InfinitePeriodicMPS(
        f, elt::Type{<:Number}, ds::AbstractVector{Int}, Ds::AbstractVector{Int}, kwargs...
    )
    return InfinitePeriodicMPS(f, elt, ComplexSpace.(ds), ComplexSpace.(Ds); kwargs...)
end

function InfinitePeriodicMPS(A::AbstractVector{<:GenericMPSTensor}; kwargs...)
    # check spaces
    leftvspaces = circshift(_firstspace.(A), -1)
    rightvspaces = conj.(_lastspace.(A))
    isnothing(findfirst(leftvspaces .!= rightvspaces)) ||
        throw(SpaceMismatch("incompatible virtual spaces $leftvspaces and $rightvspaces"))

    # check rank
    A_copy = PeriodicArray(copy.(A)) # copy to avoid side effects
    all(isfullrank, A_copy) ||
        @warn "Constructing an MPS from tensors that are not full rank"
    MPSKit.makefullrank!(A_copy)

    AR = A_copy

    leftvspaces = circshift(_firstspace.(AR), -1)
    rightvspaces = conj.(_lastspace.(AR))
    isnothing(findfirst(leftvspaces .!= rightvspaces)) ||
        throw(SpaceMismatch("incompatible virtual spaces $leftvspaces and $rightvspaces"))

    # initial guess for the gauge tensor
    V = _firstspace(A_copy[end+1])
    C₀ = isomorphism(storagetype(eltype(A_copy)), V, V)

    # initialize tensor storage
    AL = similar.(AR)
    AC = similar.(AR)
    C = similar(AR, typeof(C₀))
    ψ = InfinitePeriodicMPS{eltype(AL), eltype(C)}(AL, AR, C, AC)

    # gaugefix the MPS
    gaugefix!(ψ, A_copy, C₀; kwargs...)
    mul!.(ψ.AC, ψ.AL, ψ.C)

    return ψ
end

function InfinitePeriodicMPS(AL::AbstractVector{<:GenericMPSTensor}, C₀::MPSBondTensor; kwargs...)
    AL = PeriodicArray(copy.(AL))

    all(isfullrank, AL) ||
        @warn "Constructing an MPS from tensors that are not full rank"

    # initialize tensor storage
    AC = similar.(AL)
    AR = similar.(AL)
    T = TensorOperations.promote_contract(scalartype(AL), scalartype(C₀))
    TC = TensorOperations.tensoradd_type(T, C₀, ((1,), (2,)), false)
    C = similar(AR, TC)
    ψ = InfinitePeriodicMPS{eltype(AL), TC}(AL, AR, C, AC)

    # gaugefix the MPS
    gaugefix!(ψ, AL, C₀; order = :R, kwargs...)
    mul!.(ψ.AC, ψ.AL, ψ.C)

    return ψ
end

function InfinitePeriodicMPS(C₀::MPSBondTensor, AR::AbstractVector{<:GenericMPSTensor}; kwargs...)
    AR = PeriodicArray(copy.(AR))

    # initialize tensor storage
    AC = similar.(AR)
    AL = similar.(AR)
    T = TensorOperations.promote_contract(eltype(AR), eltype(C₀))
    TC = TensorOperations.tensoradd_type(T, C₀, ((1,), (2,)), false)
    C = similar(AR, TC)
    ψ = InfinitePeriodicMPS{eltype(AL), TC}(AL, AR, C, AC)

    # gaugefix the MPS
    gaugefix!(ψ, AR, C₀; order = :L, kwargs...)
    mul!.(ψ.AC, ψ.AL, ψ.C)

    return ψ
end

function InfinitePeriodicMPS(ψ::MPSKit.InfiniteMPS)
    return InfinitePeriodicMPS(PeriodicVector(ψ.AL), PeriodicVector(ψ.AR), PeriodicVector(ψ.C), PeriodicVector(ψ.AC))
end

#===========================================================================================
Utility
===========================================================================================#

function MPSKit.AC2(ψ::InfinitePeriodicMPS, i::Integer; kind = :ACAR)
    if kind == :ACAR
        return ψ.AC[i] * _transpose_tail(ψ.AR[i + 1])
    elseif kind == :ALAC
        return ψ.AL[i] * _transpose_tail(ψ.AC[i + 1])
    else
        throw(ArgumentError("Invalid kind: $kind"))
    end
end

Base.size(ψ::InfinitePeriodicMPS, args...) = size(ψ.AL, args...)
Base.length(ψ::InfinitePeriodicMPS) = length(ψ.AL)
Base.eltype(ψ::InfinitePeriodicMPS) = eltype(typeof(ψ))
Base.eltype(::Type{<:InfinitePeriodicMPS{A}}) where {A} = A
Base.isfinite(::Type{<:InfinitePeriodicMPS}) = false
MPSKit.GeometryStyle(::Type{<:InfinitePeriodicMPS}) = InfiniteChainStyle()

Base.copy(ψ::InfinitePeriodicMPS) = InfinitePeriodicMPS(copy(ψ.AL), copy(ψ.AR), copy(ψ.C), copy(ψ.AC))
function Base.copy!(ψ::InfinitePeriodicMPS, ϕ::InfinitePeriodicMPS)
    ψ.AL .= MPSKit._copy!!.(ψ.AL, ϕ.AL)
    ψ.AR .= MPSKit._copy!!.(ψ.AR, ϕ.AR)
    ψ.AC .= MPSKit._copy!!.(ψ.AC, ϕ.AC)
    ψ.C .= MPSKit._copy!!.(ψ.C, ϕ.C)
    return ψ
end

function Base.complex(ψ::InfinitePeriodicMPS)
    scalartype(ψ) <: Complex && return ψ
    return InfinitePeriodicMPS(complex.(ψ.AL), complex.(ψ.AR), complex.(ψ.C), complex.(ψ.AC))
end

function Base.repeat(ψ::InfinitePeriodicMPS, i::Int)
    return InfinitePeriodicMPS(repeat(ψ.AL, i), repeat(ψ.AR, i), repeat(ψ.C, i), repeat(ψ.AC, i))
end
function Base.similar(ψ::InfinitePeriodicMPS{A, B}) where {A, B}
    return InfinitePeriodicMPS{A, B}(similar(ψ.AL), similar(ψ.AR), similar(ψ.C), similar(ψ.AC))
end
function Base.circshift(ψ::InfinitePeriodicMPS, n)
    return InfinitePeriodicMPS(
        circshift(ψ.AL, n), circshift(ψ.AR, n), circshift(ψ.C, n), circshift(ψ.AC, n)
    )
end

Base.eachindex(ψ::InfinitePeriodicMPS) = eachindex(ψ.AL)
Base.eachindex(l::IndexStyle, ψ::InfinitePeriodicMPS) = eachindex(l, ψ.AL)
MPSKit.eachsite(ψ::InfinitePeriodicMPS) = PeriodicArray(eachindex(ψ))

Base.checkbounds(::Type{Bool}, ψ::InfinitePeriodicMPS, i::Integer) = true

MPSKit.site_type(::Type{<:InfinitePeriodicMPS{A}}) where {A} = A
MPSKit.bond_type(::Type{<:InfinitePeriodicMPS{<:Any, B}}) where {B} = B

MPSKit.left_virtualspace(ψ::InfinitePeriodicMPS, n::Integer) = left_virtualspace(ψ.AL[n])
MPSKit.right_virtualspace(ψ::InfinitePeriodicMPS, n::Integer) = right_virtualspace(ψ.AL[n])
MPSKit.physicalspace(ψ::InfinitePeriodicMPS, n::Integer) = physicalspace(ψ.AL[n])
MPSKit.physicalspace(ψ::InfinitePeriodicMPS) = PeriodicVector([physicalspace(ψ,n) for n in eachsite(ψ)],ψ.AL.map, ψ.AL.imap)
MPSKit.left_virtualspace(ψ::InfinitePeriodicMPS) = PeriodicVector([left_virtualspace(ψ, n) for n in eachsite(ψ)], ψ.AL.map, ψ.AL.imap)
MPSKit.right_virtualspace(ψ::InfinitePeriodicMPS) = PeriodicVector([right_virtualspace(ψ, n) for n in eachsite(ψ)], ψ.AL.map, ψ.AL.imap)

# TensorKit.space(ψ::InfinitePeriodicMPS{<:MPSTensor}, n::Integer) = space(ψ.AC[n], 2)
# function TensorKit.space(ψ::InfinitePeriodicMPS{<:GenericMPSTensor}, n::Integer)
#     t = ψ.AC[n]
#     S = spacetype(t)
#     return ProductSpace{S}(space.(Ref(t), Base.front(Base.tail(TensorKit.allind(t)))))
# end

TensorKit.norm(ψ::InfinitePeriodicMPS) = norm(ψ.AC[1])
function TensorKit.normalize!(ψ::InfinitePeriodicMPS)
    normalize!.(ψ.C)
    normalize!.(ψ.AC)
    return ψ
end

function TensorKit.dot(ψ₁::InfinitePeriodicMPS, ψ₂::InfinitePeriodicMPS; krylovdim = 30)
    # init = similar(ψ₁.AL[1], _firstspace(ψ₂.AL[1]) ← _firstspace(ψ₁.AL[1]))
    # @show init
    # randomize!(init)
    # val, = fixedpoint(
    #     flip(TransferMatrix(ψ₂.AL, ψ₁.AL)), init, :LM, MPSKit.Arnoldi(; krylovdim = krylovdim)
    # )
    init = similar(ψ₁.AL[end+1], _firstspace(ψ₂.AL[end+1]) ← _firstspace(ψ₁.AL[end+1]))
    randomize!(init)
    val, = fixedpoint(
        TransferMatrix(ψ₂.AL, ψ₁.AL), init, :LM, MPSKit.Arnoldi(; krylovdim = krylovdim)
    )
    return val
end
function Base.isapprox(ψ₁::InfinitePeriodicMPS, ψ₂::InfinitePeriodicMPS; kwargs...)
    return isapprox(dot(ψ₁, ψ₂), 1; kwargs...)
end

#===========================================================================================
Fixedpoints
===========================================================================================#

MPSKit.l_RR(ψ::InfinitePeriodicMPS, loc::Int = length(ψ) + 1) = adjoint(ψ.C[loc - 1]) * ψ.C[loc - 1]
MPSKit.l_RL(ψ::InfinitePeriodicMPS, loc::Int = 1) = ψ.C[loc - 1]
MPSKit.l_LR(ψ::InfinitePeriodicMPS, loc::Int = 1) = ψ.C[loc - 1]'

function MPSKit.l_LL(ψ::InfinitePeriodicMPS{A}, loc::Int = 1) where {A}
    return isomorphism(storagetype(A), left_virtualspace(ψ, loc), left_virtualspace(ψ, loc))
end

function MPSKit.r_RR(ψ::InfinitePeriodicMPS{A}, loc::Int = length(ψ)) where {A}
    return isomorphism(
        storagetype(A), right_virtualspace(ψ, loc), right_virtualspace(ψ, loc)
    )
end

MPSKit.r_RL(ψ::InfinitePeriodicMPS, loc::Int = length(ψ)) = ψ.C[loc]'
MPSKit.r_LR(ψ::InfinitePeriodicMPS, loc::Int = length(ψ)) = ψ.C[loc]
MPSKit.r_LL(ψ::InfinitePeriodicMPS, loc::Int = 0) = ψ.C[loc] * adjoint(ψ.C[loc])

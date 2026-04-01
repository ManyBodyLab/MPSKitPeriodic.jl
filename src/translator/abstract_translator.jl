abstract type AbstractTranslator end

# Has to be overloaded
function translate_sector end

function (tr::AbstractTranslator)(f::TensorKit.FusionTree, shift::Integer)
    new_uncoupled = map(c -> translate_sector(c, tr, shift), f.uncoupled)
    new_coupled   = translate_sector(f.coupled, tr, shift)
    return FusionTree(new_uncoupled, new_coupled, f.isdual, f.innerlines, f.vertices)
end
function (tr::AbstractTranslator)(space::TensorKit.ElementarySpace, shift::Integer)
    I = sectortype(space)
    dual = isdual(space)
    return typeof(space)(translate_sector(c, tr, shift,dual) => dim(space, c) for c in sectors(space); dual=dual)
    return Vect[I](translate_sector(c, tr, shift,dual) => dim(space, c) for c in sectors(space); dual=dual)
end
function (tr::AbstractTranslator)(space::TensorKit.ProductSpace, shift::Integer)
    isempty(space.spaces) && return space
    return ProductSpace(map(V -> tr(V, shift), space.spaces))
end
function (tr::AbstractTranslator)(space::SumSpace, shift::Int)
    return SumSpace(map(V -> tr(V, shift), space.spaces); dual = space.dual)
end
function (tr::AbstractTranslator)(space::TensorKit.HomSpace, shift::Integer)
    return TensorKit.HomSpace(tr(space.codomain, shift), tr(space.domain, shift))
end

function (tr::AbstractTranslator)(T::Union{TensorMap, BraidingTensor}, shift::Int)
    iszero(shift) && return T
    new_cod = tr(codomain(T), shift)
    new_dom = tr(domain(T), shift)
    T_new = zeros(eltype(T), new_cod, new_dom)
    for (f1, f2) in fusiontrees(T)
        f1_new = tr(f1, shift)
        f2_new = tr(f2, shift)
        T_new[f1_new, f2_new] .= T[f1, f2]
    end
    return T_new
end
function (tr::AbstractTranslator)(T::BlockTensorKit.SparseBlockTensorMap, shift::Int)
    space = tr(T.space, shift)
    data = Dict((k, tr(v, shift)) for (k, v) in T.data)
    return BlockTensorKit.SparseBlockTensorMap(data, space)
end
function (tr::AbstractTranslator)(T::BlockTensorKit.BlockTensorMap, shift::Int)
    space = tr(T.space, shift)
    data = tr.(T.data, shift)
    return BlockTensorKit.BlockTensorMap(data, space)
end


function (tr::AbstractTranslator)(T::MPSKit.JordanMPOTensor, shift::Int)
    iszero(shift) && return T
    V = tr(T.V, shift)
    A = tr(T.A, shift)
    B = tr(T.B, shift)
    C = tr(T.C, shift)
    D = tr(T.D, shift)
    return MPSKit.JordanMPOTensor(V, A, B, C, D)
end

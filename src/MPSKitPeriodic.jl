"""
Extension for MPSKit.jl that adds infinite MPS/MPO types with non-trivial boundary conditions.
"""
module MPSKitPeriodic

export InfinitePeriodicMPS
export InfinitePeriodicMPO
export InfinitePeriodicMPOHamiltonian
export FQHTranslator

import MPSKit 
using MPSKit: _transpose_tail, AbstractMPS
import MPSKit: AC2, eachsite
import MPSKit: gaugefix!
import MPSKit: GeometryStyle, InfiniteChainStyle 
import MPSKit: l_RR, l_RL, l_LR, l_LL, r_RR, r_RL, r_LR, r_LL
import MPSKit: left_virtualspace, right_virtualspace, physicalspace, bond_type, site_type
import MPSKit: GenericMPSTensor, isfullrank, _transpose_front, _firstspace, randomize!
import MPSKit: fixedpoint
import MPSKit: MPSBondTensor, AbstractTransferMatrix, MPO, MPOTensor, JordanMPOTensor
import MPSKit: MPOHamiltonian, AbstractMPSEnvironments, TransferMatrix, isidentitylevel

using BlockTensorKit: TensorMapSumSpace
using BlockTensorKit
using TensorKit
using Base.Threads: @spawn, @sync
import TensorKit: norm, normalize!, dot

using TensorOperations

import PeriodicArrays: PeriodicVector, PeriodicArray

include("states/abstractmps.jl")
include("states/infiniteperiodicmps.jl")
include("transfermatrix/transfermatrix.jl")
include("operators/infiniteperiodicmpo.jl")
include("operators/infiniteperiodicmpohamiltonian.jl")
include("environments/infinite_envs.jl")
include("environments/abstract_envs.jl")
include("algorithms/expval.jl")
include("states/ortho.jl")
include("algorithms/toolbox.jl")

include("translator/abstract_translator.jl")
include("translator/fqh_translator.jl")

end

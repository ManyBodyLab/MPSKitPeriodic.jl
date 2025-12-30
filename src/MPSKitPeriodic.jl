"""
Extension for MPSKit.jl that adds infinite MPS/MPO types with non-trivial boundary conditions.
"""
module MPSKitPeriodic

export InfinitePeriodicMPS
export InfinitePeriodicMPO
export InfinitePeriodicMPOHamiltonian

import MPSKit 
using MPSKit: _transpose_tail, AbstractMPS
import MPSKit: AC2, eachsite
import MPSKit: GeometryStyle, InfiniteChainStyle 
import MPSKit: l_RR, l_RL, l_LR, l_LL, r_RR, r_RL, r_LR, r_LL
import MPSKit: left_virtualspace, right_virtualspace, physicalspace, bond_type, site_type
import MPSKit: GenericMPSTensor, isfullrank, _left_orth!, _right_orth!, _transpose_front
import MPSKit: MPSBondTensor

using TensorKit 
import TensorKit: norm, normalize!, dot

using TensorOperations

import PeriodicArrays: PeriodicVector, PeriodicArray

include("abstractmps.jl")
include("infiniteperiodicmps.jl")
include("periodicmpo.jl")


end

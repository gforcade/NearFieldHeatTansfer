__precompile__()
module OttawaSlabCell_v3_InGaAs_InP
using FilmDataStructures, ResponseModels,eps_gold,eps_InGaAs,eps_InP
export uOttawaSlabs_v3_InGaAs_InP
# Code for initializing uOttawa slab structure 
include("uOttawaSlabStruct_v3_InGaAs_InP.jl")
precompile(uOttawaSlabs_v3_InGaAs_InP, (Float64, Float64, Int, Int, Int, Int,Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,Int))
end


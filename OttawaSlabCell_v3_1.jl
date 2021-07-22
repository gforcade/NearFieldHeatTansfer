__precompile__()
module OttawaSlabCell_v3_1
using FilmDataStructures, ResponseModels,eps_gold,eps_InAsSbP_v3,eps_InAs_ntype_v2,eps_InAs_ptype
export uOttawaSlabs_v3_1
# Code for initializing uOttawa slab structure 
include("uOttawaSlabStruct_v3_1.jl")
# Creates uOttawaSlabs(tEmit::Float64, tBck::Float64, divNCell::Int, divPCell::Int)::Tuple{lyrDsc, Array{Int64,2}} 
#precompile(uOttawaSlabs, (Float64, Float64, Int, Int))
precompile(uOttawaSlabs_v3_1, (Float64, Float64, Int, Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64,Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64))
end
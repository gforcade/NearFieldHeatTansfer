__precompile__()
module OttawaSlabCell_v2
using FilmDataStructures, ResponseModels,eps_gold,eps_InAsSbP,eps_InAs_ntype
export uOttawaSlabs_v2
# Code for initializing uOttawa slab structure 
include("uOttawaSlabStruct_v2.jl")
# Creates uOttawaSlabs(tEmit::Float64, tBck::Float64, divNCell::Int, divPCell::Int)::Tuple{lyrDsc, Array{Int64,2}} 
#precompile(uOttawaSlabs, (Float64, Float64, Int, Int))
precompile(uOttawaSlabs_v2, (Float64, Float64, Int, Int, Int, Int, Float64, Float64,Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,))
end
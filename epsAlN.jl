__precompile__()
module epsAlN
using CSV,Interpolations,DataFrames
export epsAlN_func,AlN_struct
#to model the InGaAs lattice matched to InP



struct AlNDsc
    eps1_f
    eps2_f
    T::Float64
end

@inline function AlN_struct(T)
    #get tabulated espilon data, from Adachi
    nk = CSV.read("AlN.csv",DataFrame;header=1,type=Float64 , delim="," )
    eps1 = nk[!,2]
    eps2 = nk[!,3] 
    eps_enr = nk[!,1]
    eps1_f = LinearInterpolation(eps_enr,eps1; extrapolation_bc=Flat())
    eps2_f = LinearInterpolation(eps_enr,eps2; extrapolation_bc=Flat())
    return AlNDsc(eps1_f,eps2_f,T)
end
precompile(AlN_struct,(Float64,))


@inline function  epsAlN_func(E_photon,AlNstruct)

    return  AlNstruct.eps1_f(E_photon) + im*AlNstruct.eps2_f(E_photon)  
end
precompile(epsAlN_func,(Float64,AlNDsc))


end


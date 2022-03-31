__precompile__()
module eps_InGaAs
using CSV,DataFrames,Interpolations
export InGaAs_struct, eps_InGaAs_ntype,eps_InGaAs_ptype,eps_InGaAs_imag,eVtoOmega
#to model the InGaAs lattice matched to InP
#import numpy as np
#import matplotlib.pyplot as plt

const evJ = 1.602176565e-19
const e= 1.602*10^(-19)
const eV = 1.60218*10^(-19)
const eps_0 = 8.8542*10^(-12)
const c = 3.0*10^8
const hbar = 1.05457*10^(-34) #m2kg/s
const hbar_eV = 6.5821*10^(-16) #eV s 
const m_e =9.109*10^(-31)#kg
const m_0 = 0.9109*10^(-30)
const m_0_eV = 0.51099906*10^6 #eV
const kb = 8.617333*^(10,-5) #eV K-1





#from Adachi chapter
eps_inf_InAs = 11.6
eps_inf_InSb = 15.3
eps_inf_InP = 9.9
eps_inf_GaAs = 10.86

#eV cm      #Vurgaftman , "Band parameters for III–V compound semiconductors and their alloys", 2001 
P_InAs =  9.05e-8  ##8.58e-8  ##Paskov, "Refractive indices of InSb, InAs, GaSb, InAsSb , and InGaSb: Effects of free carriers" 1997
P_GaAs =  10.47e-8
P_InSb =  9.42e-8  ##8.9e-8
P_InP  =  8.88e-8     


 #Q Parameter in Notes, No bowing parameters,  AD=InAs, BD=GaAs, AB=InGa
 @inline function  func_Q1(x,B_AD,B_BD) 
    return x*B_AD + (1.0 - x)*B_BD  
end
precompile(func_Q1,(Float64,Float64,Float64))



@inline function m_star_ptype(m_hh,m_lh)
    #calculates DOS m*_h 
    return    (m_hh^(3/2)+m_lh^(3/2))^(2/3) 
end 
precompile(m_star_ptype, (Float64,Float64))

#added this, to calculate gamma for InGaAs
@inline function Gamma_ntype_InGaAs(N_Base,T,m_star)

    #InGaAs Sotoodeh 2000
    umin = 300.0
    umax = 14000.0
    Nref = 1.3*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi = 0.48
    t1 = 1.59
    t2 = 3.68


    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    
    #this is to match with Mitapally
    if N_Base == 1.0e23
        return e/(0.041*8727.0/1e4)
    elseif N_Base == 1.0e24
        return e/(0.0483*5639.0/1e4)
    else
        return e/(m_star*mu_franc/10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
    end
    #not matching Mitapally
    #return e/(m_star*mu_franc/10000.0)

end
precompile(Gamma_ntype_InGaAs, (Float64,Float64,Float64))


@inline function Gamma_ptype_InGaAs(N_Base,T,m_star)
    
    #InGaAs: Sotoodeh 2000
    umin = 10.0
    umax = 320.0
    Nref = 4.9*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi = 0.403
    t1 = 1.59
    t2 = 0.0

    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    return e/(m_star*mu_franc/10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ptype_InGaAs, (Float64,Float64,Float64))

struct InGaAsDsc
    eps1_f
    eps2_f
    N0::Float64
    T::Float64
    E0::Float64
    gamma_ntype::Float64
    gamma_ptype::Float64
    mstar_ntype::Float64
    mstar_ptype::Float64
    mstar_ptype_lh::Float64
    mstar_ptype_hh::Float64
    epsinf::Float64
    P::Float64
    F::Float64
end

@inline function InGaAs_struct(N0,T)
    #variables for the ternary
    E0 =  0.74 #[eV] # Meghan PDF

    #get tabulated nk data, from fraunhoffer + Adachi (k: wavelength < 1000 nm)
    nk = CSV.read("InGaAs_Adachi_v2.csv",DataFrame;header=6,type=Float64 , delim="," )
    eps = (nk[!,2] + 1im*nk[!,3]).^2.0
    eps1 = reverse(real(eps))
    eps2 = reverse(imag(eps)) 
    eps_enr = reverse(1239.8 ./ nk[!,1])

    eps1_f = LinearInterpolation(eps_enr,eps1; extrapolation_bc=Flat())
    eps2_f = LinearInterpolation(eps_enr,eps2; extrapolation_bc=Flat())
    

    mstar_ntype = 0.041*m_0 #effective masses # Meghan PDF
    mstar_ptype_lh = func_Q1(0.53,0.026,0.083)*m_0
	mstar_ptype_hh = func_Q1(0.53,0.36,0.55)*m_0
    mstar_ptype = m_star_ptype(mstar_ptype_hh,mstar_ptype_lh)
    gamma_ntype = Gamma_ntype_InGaAs(N0,T,mstar_ntype)
    gamma_ptype = Gamma_ptype_InGaAs(N0,T,mstar_ptype)
    epsinf = func_Q1(0.53,eps_inf_InAs,eps_inf_GaAs)
    P = func_Q1(0.53,P_InAs,P_GaAs)
    Nc = 2.0*(3.0*E0*kb*T/(8.0*pi*P^2.0))^(3.0/2.0) 
    #calculating the fermi energy
	#res = optimize(t -> fermi(t,E0/(kb*T),pi^(1/2)/2*N0/(Nc*10.0^6)),[E0],Newton()) #Nc from cm-3 to m-3
	#F = Optim.minimizer(res)[1]*kb*T + E0 #fermi energy relative to the valence band
    F = 0.0
    return InGaAsDsc(eps1_f,eps2_f,N0,T,E0,gamma_ntype ,gamma_ptype ,mstar_ntype ,mstar_ptype,mstar_ptype_lh,mstar_ptype_hh,epsinf,P,F)
end
precompile(InGaAs_struct,(Float64,Float64,))

@inline function eVtoOmega(energy) #energy in eV
    return(energy*evJ/hbar)
end
precompile(eVtoOmega, (Float64,))



@inline function  eps_InGaAs_ntype(E_photon,InGaAsstruct)

    return  InGaAsstruct.eps1_f(E_photon) + im*InGaAsstruct.eps2_f(E_photon)  + epsFCL_InGaAs(eVtoOmega(E_photon),InGaAsstruct.mstar_ntype,InGaAsstruct.gamma_ntype,InGaAsstruct.N0,InGaAsstruct.epsinf)
end
precompile(eps_InGaAs_ntype,(Float64,InGaAsDsc))


@inline function  eps_InGaAs_ptype(E_photon,InGaAsstruct)
    
    return  InGaAsstruct.eps1_f(E_photon) + im*InGaAsstruct.eps2_f(E_photon)  + epsFCL_InGaAs(eVtoOmega(E_photon),InGaAsstruct.mstar_ptype,InGaAsstruct.gamma_ptype,InGaAsstruct.N0,InGaAsstruct.epsinf)

end
precompile(eps_InGaAs_ptype,(Float64,InGaAsDsc))



@inline function eps_InGaAs_imag(enr,InGaAsstruct)


    if(enr>InGaAsstruct.E0)
        
        return  InGaAsstruct.eps2_f(enr)
    else
        return 0.0
    end
end
precompile(eps_InGaAs_imag,(Float64,InGaAsDsc))








@inline function epsFCL_InGaAs(omega,mstar,gamma,N_base,eps_inf)
    #calculating the free carrier and lattice effects for InP

    #eps,gamma parameters from Sadao Adachi: Optical Properties of Crystalline and Amorphous Semiconductors (1999)
    #To and LO parameters: Adachi - Properties of semiconductor alloys
    #GaAs
    g_GaAs = 4.52*10.0^11  #s^(-1), gamma
    o_TO_GaAs = 5.05*10.0^(13) #s^(-1)
    o_LO_GaAs = 5.50*10.0^(13)

    #InAs
    #eps_inf_InAs = 12.3
    g_InAs = 9.23*10.0^11  #s^(-1), gamma 
    o_TO_InAs = 4.14*10.0^(13) #s^(-1) 
    o_LO_InAs = 4.55*10.0^(13)
    
    o_p_square=N_base*e^2/(eps_0*eps_inf*mstar) #s^(-2)
    
    #first term is free carrier effect, the rest are the lattice effects, weiighted by the atom fraction
    return eps_inf*(- o_p_square/(omega*(omega+ im*gamma)) + 0.53*(o_LO_InAs^2-o_TO_InAs^2)/(o_TO_InAs^2-omega^2-im*omega*g_InAs) + 0.47*(o_LO_GaAs^2-o_TO_GaAs^2)/(o_TO_GaAs^2-omega^2-im*omega*g_GaAs))
end
precompile(epsFCL_InGaAs,(Float64,Float64,Float64,Float64))


end


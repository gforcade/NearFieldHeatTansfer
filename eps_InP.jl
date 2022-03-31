__precompile__()
module eps_InP
using CSV,DataFrames,Interpolations
export InP_struct, eps_InP_ntype,eps_InP_ptype,eps_InP_imag,eVtoOmega
#material properties for InP

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





#below are the new important parameters, everything else is not used
eps_inf_InAs = 11.6
eps_inf_InSb = 15.3
eps_inf_InP = 9.9

#eV cm      #Vurgaftman , "Band parameters for III–V compound semiconductors and their alloys", 2001 
P_InAs =  9.05e-8  ##8.58e-8  ##Paskov, "Refractive indices of InSb, InAs, GaSb, InAsSb , and InGaSb: Effects of free carriers" 1997
P_InSb =  9.42e-8  ##8.9e-8
P_InP  =  8.88e-8     






@inline function m_star_ptype(m_hh,m_lh)
    #calculates DOS m*_h 
    return    (m_hh^(3/2)+m_lh^(3/2))^(2/3) 
end 
precompile(m_star_ptype, (Float64,Float64))

#added this, to calculate gamma for InP
@inline function Gamma_ntype_InP(N_Base,T,m_star)

    #InP: Sotoodeh 2000
    umin = 400.0
    umax = 5200.0
    Nref = 3.0*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi = 0.47
    t1 = 2.0
    t2 = 3.25


    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    
    #this is to match with Mitapally
    if N_Base == 1.0e24
        return e/(0.0831*2096.0/1e4)
    else
        return e/(m_star*mu_franc/10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
    end

    #not matching Mitapally
    #return e/(m_star*mu_franc/10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ntype_InP, (Float64,Float64,Float64))


@inline function Gamma_ptype_InP(N_Base,T,m_star)
    
    """
    #not matching Mitapally
    #InP: Sotoodeh 2000
    umin = 10.0
    umax = 170.0
    Nref = 4.97*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi = 0.62
    t1 = 2.0
    t2 = 2.0
    """
    
    #InP: matching Mitapally 2021
    umin = 0.0
    umax = 150.0
    Nref = 2.0*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi = 0.5
    t1 = 2.0
    t2 = 2.0
    

    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    return e/(m_star*mu_franc/10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ptype_InP, (Float64,Float64,Float64))

struct InPDsc
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

@inline function InP_struct(N0,T)
    E0 =  1.344 #[eV] # Meghan PDF
    #get n k values from tabulated data
    nk = CSV.read("InP_nk_ioffe.par",DataFrame;header=6,type=Float64 , delim="\t" )
    eps = (nk[!,2] + 1im*nk[!,3]).^2.0
    eps1 = reverse(real(eps))
    eps2 = reverse(imag(eps)) 
    eps_enr = reverse(1239.8 ./ nk[!,1])

    eps1_f = LinearInterpolation(eps_enr,eps1; extrapolation_bc=Flat())
    eps2_f = LinearInterpolation(eps_enr,eps2; extrapolation_bc=Flat())


    mstar_ntype = 0.079*m_0 #effective masses # Meghan PDF
    mstar_ptype_lh = 0.11*m_0 # Meghan PDF
	mstar_ptype_hh = 0.69*m_0 # Meghan PDF
    #mstar_ptype = m_star_ptype(mstar_ptype_hh,mstar_ptype_lh)  ####not matching Mitapally
    mstar_ptype = 0.6*m_0 #from Mitapally
    gamma_ntype = Gamma_ntype_InP(N0,T,mstar_ntype)
    gamma_ptype = Gamma_ptype_InP(N0,T,mstar_ptype)
    epsinf = eps_inf_InP
    P = P_InP
    F = 0.0
    return InPDsc(eps1_f,eps2_f,N0,T,E0,gamma_ntype ,gamma_ptype ,mstar_ntype ,mstar_ptype,mstar_ptype_lh,mstar_ptype_hh,epsinf,P,F)
end
precompile(InP_struct,(Float64,Float64,))

@inline function eVtoOmega(energy) #energy in eV
    return(energy*evJ/hbar)
end
precompile(eVtoOmega, (Float64,))



@inline function  eps_InP_ntype(E_photon,InPstruct)

    return  InPstruct.eps1_f(E_photon) + im*InPstruct.eps2_f(E_photon)  + epsFCL_InP(eVtoOmega(E_photon),InPstruct.mstar_ntype,InPstruct.gamma_ntype,InPstruct.N0,InPstruct.epsinf)
end
precompile(eps_InP_ntype,(Float64,InPDsc))


@inline function  eps_InP_ptype(E_photon,InPstruct)
    
    return  InPstruct.eps1_f(E_photon) + im*InPstruct.eps2_f(E_photon)  + epsFCL_InP(eVtoOmega(E_photon),InPstruct.mstar_ptype,InPstruct.gamma_ptype,InPstruct.N0,InPstruct.epsinf)

end
precompile(eps_InP_ptype,(Float64,InPDsc))



@inline function eps_InP_imag(enr,InPstruct)


    if(enr>InPstruct.E0)

        return  InPstruct.eps2_f(enr)
    else
        return 0.0
    end
end
precompile(eps_InP_imag,(Float64,InPDsc))








@inline function epsFCL_InP(omega,mstar,gamma,N_base,eps_inf)
    #calculating the free carrier and lattice effects for InP

    #eps,gamma parameters from Sadao Adachi: Optical Properties of Crystalline and Amorphous Semiconductors (1999)
    #To and LO parameters: Adachi - Properties of semiconductor alloys
    #InP
    g_InP = 3.58*10.0^11  #s^(-1), gamma
    o_TO_InP = 5.74*10.0^(13) #s^(-1)
    o_LO_InP = 6.52*10.0^(13)
    
    o_p_square=N_base*e^2/(eps_0*eps_inf*mstar) #s^(-2)
    
    #first term is free carrier effect, the rest are the lattice effects, weiighted by the atom fraction
    return eps_inf*(- o_p_square/(omega*(omega+ im*gamma))  + (o_LO_InP^2-o_TO_InP^2)/(o_TO_InP^2-omega^2-im*omega*g_InP))
end
precompile(epsFCL_InP,(Float64,Float64,Float64,Float64))
end
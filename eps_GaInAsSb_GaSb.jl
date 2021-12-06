__precompile__()
module eps_InAsSbP_v3
using eps_InAs_ntype_v2,Optim
export G_E1_E2,func_Q,func_Q2,func_f,InAsSbP_struct, eps_InAsSbP_xy_ntype,eps_InAsSbP_xy_ptype,eps_InAsSbP_imag_xy_ntype,eps_InAsSbP_imag_xy_ptype,Gamma_ntype_InAsSbP,Gamma_ptype_InAsSbP,eVtoOmega,epsFCL_InAsSbP
#Used this version in the end, accurate imaginary part from the range of 0.1 to 1.5, above 1.5, imaginary part becomes negative"""
#tried different model dielectric function. Don't use this.


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
const T_global = 300.0
const kb = 8.617333*^(10,-5) #eV K-1


#below are the new important parameters, everything else is not used
eps_inf_GaAs = 10.86
eps_inf_GaSb = 14.2
eps_inf_InAs = 11.6
eps_inf_InSb = 15.3


#Vurgaftman , "Band parameters for III–V compound semiconductors and their alloys", 2001 
P_GaAs = 10.47e-8   #[eV cm]      
P_GaSb =  10.14e-8
P_InAs =  9.05e-8  
P_InSb =  9.42e-8  
P_InP  =  8.88e-8     




#Q Parameter in Notes
@inline function  func_Q(x,y,B_AD,B_BD,B_AC,B_CD,C_ABC)
    return x*B_AD + y*B_BD + (1-x-y)*B_CD  + C_ABC*x*y*(1-x-y)
end
precompile(func_Q,(Float64,Float64,Float64,Float64,Float64,Float64,Float64))

 #Ga(x)In(1-x)As(y)Sb(1-y) No bowing parameters,  AC=GaAs, AD=GaSb, BC=InAs BD=InSb
@inline function  func_Q2(x,y,B_AC,B_AD,B_BC,B_BD) 
    return x*y*B_AC + x*(1-y)*B_AD + (1-x)*y*B_BC + (1-x)*(1-y)*B_BD 
end
precompile(func_Q2,(Float64,Float64,Float64,Float64,Float64,Float64))


 #Ga(x)In(1-x)As(y)Sb(1-y) No bowing parameters, Mathiessens rule ouput value (not 1/value),  AC=GaAs, AD=GaSb, BC=InAs BD=InSb
 @inline function  func_Q3(x,y,B_AC,B_AD,B_BC,B_BD) 
    return 1/(x*y/B_AC + x*(1-y)/B_AD + (1-x)*y/B_BC + (1-x)*(1-y)/B_BD)
end
precompile(func_Q3,(Float64,Float64,Float64,Float64,Float64,Float64))


@inline function  func_f(z)
    return z^(-2)*(2.0-(1.0+z)^0.5-(1.0-z)^0.5)
end
precompile(func_f,(Float64,))

@inline function m_star_ptype(m_hh,m_lh)
    #calculates DOS m*_h 
    return    (m_hh^(3/2)+m_lh^(3/2))^(2/3) 
end 
precompile(m_star_ptype, (Float64,Float64))

#added this, to calculate gamma for InAsSbP
@inline function Gamma_ntype_GaInAsSb(x,y,N_Base,T,m_star)
    #mobility params from "Numerical analysis of the short-circuit current density in GaInAsSb thermophotovoltaic diodes" Peng et al. 2009
    #parameter t from "Simulation of temperature-dependent material parameters and device performances for GaInAsSb thermophotovoltaic cell" Peng et al. 2011


    #GaAs
    umin_GaAs = 500.0
    umax_GaAs = 9400.0
    Nref_GaAs = 0.6*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_GaAs = 0.3940
    t1_GaAs = 2.1               #from
    t2_GaAs = 3.0

    #GaSb
    umin_GaSb = 1050.0
    umax_GaSb = 5650.0
    Nref_GaSb = 2.8*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_GaSb = 1.05
    t1_GaSb = 2.0
    t2_GaSb = 2.8

    #InAs 
    umin_InAs = 1000.0
    umax_InAs = 34000.0
    Nref_InAs = 1.1*10^18*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InAs = 0.32
    t1_InAs = 1.57
    t2_InAs = 3.0

    #InSb
    umin_InSb = 5000.0
    umax_InSb = 78000.0
    Nref_InSb = 0.7*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InSb = 0.7
    t1_InSb = 2.0
    t2_InSb = 2.8


    #interpolated values
    umin = func_Q3(x,y,umin_GaAs,umin_GaSb,umin_InAs,umin_InSb)
    umax = func_Q3(x,y,umax_GaAs,umax_GaSb,umax_InAs,umax_InSb)
    Nref = func_Q2(x,y,Nref_GaAs,Nref_GaSb,Nref_InAs,Nref_InSb)
    phi = func_Q2(x,y,phi_GaAs,phi_GaSb,phi_InAs,phi_InSb)
    t1 = func_Q2(x,y,t1_GaAs,t1_GaSb,t1_InAs,t1_InSb)
    t2 = func_Q2(x,y,t2_GaAs,t2_GaSb,t2_InAs,t2_InSb)

    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)    #from Francoeur 
    return e/(m_star*mu_franc*10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ntype_GaInAsSb, (Float64,Float64,Float64,Float64,Float64))


@inline function Gamma_ptype_GaInAsSb(x,y,N_Base,T,m_star)
     #mobility params from "Numerical analysis of the short-circuit current density in GaInAsSb thermophotovoltaic diodes" Peng et al. 2009
    #parameter t from "Simulation of temperature-dependent material parameters and device performances for GaInAsSb thermophotovoltaic cell" Peng et al. 2011


    #GaAs
    umin_GaAs = 20.0
    umax_GaAs = 491.5
    Nref_GaAs = 1.48*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_GaAs = 0.38
    t1_GaAs = 2.2               #from
    t2_GaAs = 3.0

    #GaSb
    umin_GaSb = 190.0
    umax_GaSb = 875.0
    Nref_GaSb = 9.0*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_GaSb = 0.65
    t1_GaSb = 1.7
    t2_GaSb = 2.7

    #InAs 
    umin_InAs = 20.0
    umax_InAs = 530.0
    Nref_InAs = 1.1*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InAs = 0.46
    t1_InAs = 2.3
    t2_InAs = 3.0

    #InSb
    umin_InSb = 100.0
    umax_InSb = 750.0
    Nref_InSb = 6.0*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InSb = 0.6
    t1_InSb = 1.7
    t2_InSb = 2.7


    #interpolated values
    umin = func_Q3(x,y,umin_GaAs,umin_GaSb,umin_InAs,umin_InSb)
    umax = func_Q3(x,y,umax_GaAs,umax_GaSb,umax_InAs,umax_InSb)
    Nref = func_Q2(x,y,Nref_GaAs,Nref_GaSb,Nref_InAs,Nref_InSb)
    phi = func_Q2(x,y,phi_GaAs,phi_GaSb,phi_InAs,phi_InSb)
    t1 = func_Q2(x,y,t1_GaAs,t1_GaSb,t1_InAs,t1_InSb)
    t2 = func_Q2(x,y,t2_GaAs,t2_GaSb,t2_InAs,t2_InSb)

    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    return e/(m_star*mu_franc*10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ptype_GaInAsSb, (Float64,Float64,Float64,Float64,Float64))

struct GaInAsSbDsc
    x::Float64
    y::Float64
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

@inline function GaInAsSb_struct(x,y,N0,T)
    #provides various parameters for given material composition, temperature and doping

    E0 =  0.28 - 0.16*x + 0.60*x^2  #only good for GaInAsSb lattice matched to GaSb #from Adachi 2017    
    mstar_ntype = func_Q2(x,y,0.067,0.039,0.024,0.013)*m_0 #effective masses from Adachi, "Properties of Semiconductor alloys: ...," 2009
    mstar_ptype_lh = func_Q2(x,y,0.083,0.043,0.026,0.014)*m_0
	mstar_ptype_hh = func_Q2(x,y,0.55,0.37,0.36,0.38)*m_0
    mstar_ptype = m_star_ptype(mstar_ptype_hh,mstar_ptype_lh)
    gamma_ntype = Gamma_ntype_GaInAsSb(x,y,N0,T,mstar_ntype)
    gamma_ptype = Gamma_ptype_GaInAsSb(x,y,N0,T,mstar_ptype)
    epsinf = func_Q2(x,y,eps_inf_GaAs,eps_inf_GaSb,eps_inf_InAs,eps_inf_InSb)
    P = func_Q2(x,y,P_GaAs,P_GaSb,P_InAs,P_InSb)
    Nc = 2.0*(3.0*E0*kb*T/(8.0*pi*P^2.0))^(3.0/2.0) 
    #calculating the fermi energy
	res = optimize(t -> fermi(t,E0/(kb*T),pi^(1/2)/2*N0/(Nc*10.0^6)),[E0],Newton()) #Nc from cm-3 to m-3
	F = Optim.minimizer(res)[1]*kb*T + E0 #fermi energy relative to the valence band
    return GaInAsSbDsc(x,y,N0,T,E0,gamma_ntype ,gamma_ptype ,mstar_ntype ,mstar_ptype,mstar_ptype_lh,mstar_ptype_hh,epsinf,P,F)
end
precompile(GaInAsSb_struct,(Float64,Float64,Float64,Float64,))

@inline function eVtoOmega(energy) #energy in eV
    return(energy*evJ/hbar)
end
precompile(eVtoOmega, (Float64,))



@inline function  eps_GaInAsSb_xy_ntype(E_photon,GaInAsSbstruct)

    #adding IB or no IB. IB - epsinf since I don't want to double count
    #if(E_photon>InAsSbPstruct.E0_InAsSbP)
    return epsIB(eVtoOmega(E_photon),GaInAsSbstruct.N0,GaInAsSbstruct.T,GaInAsSbstruct.E0,GaInAsSbstruct.epsinf,GaInAsSbstruct.P,GaInAsSbstruct.mstar_ptype_hh,GaInAsSbstruct.mstar_ptype_lh,GaInAsSbstruct.F,0.0) - GaInAsSbstruct.epsinf + epsFCL_InAsSbP(GaInAsSbstruct.x,GaInAsSbstruct.y,eVtoOmega(E_photon),GaInAsSbstruct.mstar_ntype,GaInAsSbstruct.gamma_ntype,GaInAsSbstruct.N0,GaInAsSbstruct.epsinf)
    #else
    #    return epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(E_photon),InAsSbPstruct.mstar_ntype,InAsSbPstruct.gamma_ntype,InAsSbPstruct.N0,InAsSbPstruct.epsinf)
    #end
end
precompile(eps_GaInAsSb_xy_ntype,(Float64,GaInAsSbDsc))


@inline function  eps_GaInAsSb_xy_ptype(E_photon,GaInAsSbstruct)
    
    #adding the IB and FCL to eps
    #if(E_photon>InAsSbPstruct.E0_InAsSbP)
    return epsIB(eVtoOmega(E_photon),GaInAsSbstruct.N0,GaInAsSbstruct.T,GaInAsSbstruct.E0,GaInAsSbstruct.epsinf,GaInAsSbstruct.P,GaInAsSbstruct.mstar_ptype_hh,GaInAsSbstruct.mstar_ptype_lh,GaInAsSbstruct.F,1.0) - GaInAsSbstruct.epsinf + epsFCL_GaInAsSb(GaInAsSbstruct.x,GaInAsSbstruct.y,eVtoOmega(E_photon),GaInAsSbstruct.mstar_ptype,GaInAsSbstruct.gamma_ptype,GaInAsSbstruct.N0,GaInAsSbstruct.epsinf)
    #else
    #    return epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(E_photon),InAsSbPstruct.mstar_ptype,InAsSbPstruct.gamma_ptype,InAsSbPstruct.N0,InAsSbPstruct.epsinf)
    #end
end
precompile(eps_GaInAsSb_xy_ptype,(Float64,GaInAsSbDsc))



@inline function eps_GaInAsSb_imag_xy_ptype(enr,GaInAsSbstruct)
    #In the code InAs_xSb_yP_1-x-y,  from III-V Ternary and Quaternary Compounds Table 30.14 uses InP_xAs_ySb_1-x-y, and E0 = 0.512+0.030*y-0.183*y^2
    #E_bandgap = 0.512 + 0.030*InAsSbPstruct.x-0.183*InAsSbPstruct.x^2   
    if(enr>GaInAsSbstruct.E0)
        return imag(epsIB(eVtoOmega(E_photon),GaInAsSbstruct.N0,GaInAsSbstruct.T,GaInAsSbstruct.E0,GaInAsSbstruct.epsinf,GaInAsSbstruct.P,GaInAsSbstruct.mstar_ptype_hh,GaInAsSbstruct.mstar_ptype_lh,GaInAsSbstruct.F,1.0))
    else
        return 0.0
    end
end
precompile(eps_GaInAsSb_imag_xy_ptype,(Float64,GaInAsSbDsc))




@inline function eps_GaInAsSb_imag_xy_ntype(enr,GaInAsSbstruct)
    #In the code InAs_xSb_yP_1-x-y,  from III-V Ternary and Quaternary Compounds Table 30.14 uses InP_xAs_ySb_1-x-y, and E0 = 0.512+0.030*y-0.183*y^2
    #E_bandgap = 0.512 + 0.030*InAsSbPstruct.x-0.183*InAsSbPstruct.x^2   
    if(enr>GaInAsSbstruct.E0)
        return imag(epsIB(eVtoOmega(E_photon),GaInAsSbstruct.N0,GaInAsSbstruct.T,GaInAsSbstruct.E0,GaInAsSbstruct.epsinf,GaInAsSbstruct.P,GaInAsSbstruct.mstar_ptype_hh,GaInAsSbstruct.mstar_ptype_lh,GaInAsSbstruct.F,0.0))
    else
        return 0.0
    end
end
precompile(eps_GaInAsSb_imag_xy_ntype,(Float64,GaInAsSbDsc))



@inline function epsFCL_GaInAsSb(x,y,omega,mstar,gamma,N_base,eps_inf)
    #calculating the free carrier and lattice effects for InAsSbP

    #eps,gamma parameters from Sadao Adachi: Optical Properties of Crystalline and Amorphous Semiconductors (1999)
    #To and LO parameters: Adachi - Properties of semiconductor alloys
    #InAs
    #eps_inf_InAs = 12.3
    g_InAs = 9.23*10.0^11  #s^(-1), gamma 
    o_TO_InAs = 4.14*10.0^(13) #s^(-1) 
    o_LO_InAs = 4.55*10.0^(13)


    #InSb
    #eps_inf_InSb = 15.68
    g_InSb = 5.41*10.0^11  #s^(-1), gamma
    o_TO_InSb = 3.38*10.0^(13) #s^(-1)
    o_LO_InSb =3.59*10.0^(13)


    #InP
    #eps_inf_InP = 9.66
    g_InP = 3.58*10.0^11  #s^(-1), gamma
    o_TO_InP = 5.74*10.0^(13) #s^(-1)
    o_LO_InP = 6.52*10.0^(13)
    
    
    o_p_square=N_base*e^2/(eps_0*eps_inf*mstar) #s^(-2)
    
    #first term is free carrier effect, the rest are the lattice effects, weiighted by the atom fraction
    return eps_inf*(1.0 - o_p_square/(omega*(omega+ im*gamma)) + x*(o_LO_InAs^2-o_TO_InAs^2)/(o_TO_InAs^2-omega^2-im*omega*g_InAs)+y*(o_LO_InSb^2-o_TO_InSb^2)/(o_TO_InSb^2-omega^2-im*omega*g_InSb)+(1-x-y)*(o_LO_InP^2-o_TO_InP^2)/(o_TO_InP^2-omega^2-im*omega*g_InP))
end
precompile(epsFCL_InAsSbP,(Float64,Float64,Float64,Float64,Float64,Float64))
end
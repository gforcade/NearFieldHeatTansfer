__precompile__()
module eps_InAsSbP_v3
using eps_InAs_ntype_v2,Optim
export G_E1_E2,func_Q,func_Q2,E_T,func_f,InAsSbP_struct, eps_InAsSbP_xy_ntype,eps_InAsSbP_xy_ptype,eps_InAsSbP_imag_xy_ntype,eps_InAsSbP_imag_xy_ptype,Gamma_ntype_InAsSbP,Gamma_ptype_InAsSbP,eVtoOmega,epsFCL_InAsSbP
#Used this version in the end, accurate imaginary part from the range of 0.1 to 1.5, above 1.5, imaginary part becomes negative"""
#tried different model dielectric function. Don't use this.
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
const T_global = 300.0
const kb = 8.617333*^(10,-5) #eV K-1

#list of constants: list of constants for InAsSbP, AD=InAs, BD=InSb, CD=InP ABC=AsSbP  #Cuevas Paper Table1
#can't find bowing constants for AsSbP
#Bowing Constants for GaIn, AsSb from ZhangThermoConversion
A_InAs = 19.950   #A (eV^1.5) for eps_1, Cuevas
A_InSb = 16.830   #A (eV^1.5) for eps_1, Cueva #A_InP =  6.57     #Need fitting Adachi Table1
A_InP =  4.062    # Modeling the optical constants of GaP, InP, and InAs - (1999)
C_InAs = 1.750  #C(Dimension-less) for eps_3 Cuevas
C_InSb = 2.64
C_InP=  1.49#0.5600889435237922
D_InAs =21.511 #D(Dimension-less) for eps_4 Cuevas
D_InSb = 38.830 #D_InP = 60.4 #Need fitting Adachi Table1
D_InP=  60.4#-147.6161185407525
G_E0_InAs = 8.660 #in eV Cuevas
G_E0_InSb = 4.20
G_E0_InP= 0.036 #0.12337661235349834

#InAs
epsinf_InAs = 1.142
A_InAs = 3.264
G_E0_InAs = 0.165
alpha0_InAs = 4.111
B1_InAs = 2.704
B2_InAs = 2.765
E0_InAs = 0.356
E0Delta0_InAs = 0.736
G_E1_InAs = 4.394
alpha1_InAs = 0.027
f2_InAs = 3.984
G_E2_InAs = 0.619
alpha2_InAs = 0.033



@inline function  G_E1_E2(gamma_L,g,T)
    return gamma_L + g*T
end
precompile(G_E1_E2,(Float64,Float64,Float64))


G_E1_InAs = G_E1_E2(0.25,0.3,T_global)*0.001 #gamma in meV, change to eV #Cuevas Table 3
G_E1_InSb = G_E1_E2(0.15,0.28,T_global)*0.001
G_E1_InP  = 5.217 # Adachi Table1

G_E2_InAs =  G_E1_E2(403,0.3,T_global)*0.001 #gamma in meV, change to eV #Cuevas Table 3
G_E2_InSb =  G_E1_E2(633.9,0.49,T_global)*0.001 
G_E2_InP =  0.10 # Adachi Table1

G_Eg_InAs =  0.350#(eV) Cuevas
G_Eg_InSb = 0.56 
#Cuevas Table2, The energies are in eV, NO bowing constants for E0??Only temperature depence?
E0_InAs = 0.42  #at 300K, E0_T = 1.5,at 0K 1.4
dt0_InAs = 0.28*0.001
b0_InAs = 93.0
E0_InSb = 0.24
dt0_InSb = 0.32*0.001
b0_InSb = 170.0 #E0_InP  = 1.35 #Need fitting Adachi Table1
E0_InP = 1.344  #at 300K, E0_T = 1.5,at 0K 1.4
#dt0_InP = 0.49*0.001
#b0_InP = 327.0

E0Delta0_InAs = 0.79#-0.42 
dtE0Delta0_InAs = 0.34*0.001
bE0Delta0_InAs  = 248.0
E0Delta0_InSb = 1.20#-0.24
dtE0Delta0_InSb = 0.32*0.001
bE0Delta0_InSb = 170.0 #E0Delta0_GaAs = 1.85#-0.152 dtE0Delta0_GaAs = 0.35*0.001 bE0Delta0_GaAs = 225.0 E0Delta0_GaSb = 1.57#-0.81 dtE0Delta0_GaSb = 0.42*0.001
E0Delta0_InP = 1.414 # Adachi Table B15-1



E1_InAs = 2.5 #Cuevas Table2 
dt1_InAs = 0.5*0.001
b1_InAs = 0.0
E1_InSb = 2.0 
dt1_InSb = 0.68*0.001
b1_InSb = 132.0
E1_InP = 3.17 # Adachi Table B15-1
  
E1Delta1_InAs =2.88#-2.61 #Cuevas Table 2
dtE1Delta1_InAs = 0.5*0.001
bE1Delta1_InAs = 0.0
E1Delta1_InSb = 2.49#-2.00
dtE1Delta1_InSb = 0.65*100
bE1Delta1_InSb = 170.0 
E1Delta1_InP = 3.29 #Adachi  Table B15-1
  
Delta1_InAs = 2.88-2.61 #Cuevas
Delta1_InSb = 2.49-2.00 #Delta1_InP = 0.15 #Adachi  Table1
Delta1_InP =  3.29-3.17 

E2_InAs = 4.70 #from Zhang Thermo Conversion
dt2_InAs = 0.56*0.001
b2_InAs = 0.0
E2_InSb = 3.9
dt2_InSb = 0.54*0.001
b2_InSb = 0.0
E2_InP = 5.1 #Adachi  Table B15-1

Eg_InAs = 1.13 #Cuevas Table 2
dtg_InAs = 0.28*0.001
bg_InAs = 93.0
Eg_InSb = 0.93
dtg_InSb = 0.32*0.001
bg_InSb = 170.0
Eg_InP  = 2.05 #Adachi Table B15-1

#Cuevas paper Table IV in angstrong
alpha_lc_InAs = 6.06
alpha_lc_InSb = 6.48
alpha_lc_InP =  5.87 
#alpha_lc_GaAs = 5.65 alpha_lc_GaSb = 6.10 alpha_lc_GaIn = 0.055  #Cuevas no value #zhang 0.055 alpha_lc_AsSb = -0.023 #Cuevas no value #zhang -0.023

#Fitted Parameters
G_Eg_InP=  0.8902957174564856



#below are the new important parameters, everything else is not used
eps_inf_InAs = 12.3
eps_inf_InSb = 15.68
eps_inf_InP = 9.66

#eV cm      #Vurgaftman , "Band parameters for III–V compound semiconductors and their alloys", 2001 
P_InAs =  9.05e-8  ##8.58e-8  ##Paskov, "Refractive indices of InSb, InAs, GaSb, InAsSb , and InGaSb: Effects of free carriers" 1997
P_InSb =  9.42e-8  ##8.9e-8
P_InP  =  8.88e-8     

#cm-3 #all from sentaurus
Nc_InAs = 9.3301e16 
Nc_InSb = 3.7195e16
Nc_InP  = 5.6e17


#Q Parameter in Notes
@inline function  func_Q(x,y,B_AD,B_BD,B_AC,B_CD,C_ABC)
    return x*B_AD + y*B_BD + (1-x-y)*B_CD  + C_ABC*x*y*(1-x-y)
end
precompile(func_Q,(Float64,Float64,Float64,Float64,Float64,Float64,Float64))

 #Q Parameter in Notes, No bowing parameters,  AD=InAs, BD=InSb, CD=InP ABC=AsSbP
@inline function  func_Q2(x,y,B_AD,B_BD,B_CD) 
    return x*B_AD + y*B_BD + (1-x-y)*B_CD 
end
precompile(func_Q2,(Float64,Float64,Float64,Float64,Float64))

@inline function  E_T(E0,dt,beta)
    return E0-dt*300.0^2/(300.0+beta)
end
precompile(E_T,(Float64,Float64,Float64))

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
@inline function Gamma_ntype_InAsSbP(x,y,N_Base,T,m_star)
    
    #InAs 
    #me_InAs = 0.024*m_0 #Adachi "Properties of Semiconductors alloys," 2009
    umin_InAs = 1000.0
    umax_InAs = 34000.0
    Nref_InAs = 1.1*10^18*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InAs = 0.32
    t1_InAs = 1.57
    t2_InAs = 3.0

    #InSbP (same as Sentaurus)
    #me_InSbP = 0.063*m_0 #Adachi "Properties of Semiconductors alloys," 2009
    umin_InSbP = 559.6
    umax_InSbP = 6756.9
    Nref_InSbP = 2.287*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InSbP = 0.5413
    t1_InSbP = 2.0
    t2_InSbP = 3.1105

    #InP

    #interpolated values
    #m_star = x*me_InAs + (1-x)*me_InSbP
    umin = (x/umin_InAs + (1-x)/umin_InSbP)
    umax = (x/umax_InAs + (1-x)/umax_InSbP)
    Nref = x*Nref_InAs + (1-x)*Nref_InSbP
    phi = x*phi_InAs + (1-x)*phi_InSbP
    t1 = x*t1_InAs + (1-x)*t1_InSbP
    t2 = x*t2_InAs + (1-x)*t2_InSbP

    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    return e/(m_star*mu_franc*10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ntype_InAsSbP, (Float64,Float64,Float64,Float64,Float64))


@inline function Gamma_ptype_InAsSbP(x,y,N_Base,T,m_star)
    
    #InAs
    #me_InAs = 0.41*m_0
    umin_InAs = 1000.0
    umax_InAs = 34000.0
    Nref_InAs = 1.1*10^18*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InAs = 0.32
    t1_InAs = 1.57
    t2_InAs = 3.0

    #InSbP (same as Sentaurus)
    #me_InSbP = 0.54*m_0
    umin_InSbP = 13.87
    umax_InSbP = 208.8
    Nref_InSbP = 5.2203*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InSbP = 0.6138
    t1_InSbP = 1.907
    t2_InSbP = 2.907

    #interpolated values
    #m_star = x*me_InAs + (1-x)*me_InSbP
    umin = (x/umin_InAs + (1-x)/umin_InSbP)
    umax = (x/umax_InAs + (1-x)/umax_InSbP)
    Nref = x*Nref_InAs + (1-x)*Nref_InSbP
    phi = x*phi_InAs + (1-x)*phi_InSbP
    t1 = x*t1_InAs + (1-x)*t1_InSbP
    t2 = x*t2_InAs + (1-x)*t2_InSbP

    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    return e/(m_star*mu_franc*10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ptype_InAsSbP, (Float64,Float64,Float64,Float64,Float64))

struct InAsSbPDsc
    x::Float64
    y::Float64
    N0::Float64
    T::Float64
    E0_InAsSbP::Float64
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

@inline function InAsSbP_struct(x,y,N0,T)
     #assume equal proportion of As,Sb and P
    #list of constants for InAsSbP, AD=InAs, BD=InSb, CD=InP ABC=AsSbP  #Cuevas Paper Table1
    #Bowing Constants for GaIn, AsSb from ZhangThermoConversion
    #=
	E0_T_InAs = E_T(E0_InAs,dt0_InAs,b0_InAs)
	E0_T_InSb = E_T(E0_InSb,dt0_InSb,b0_InSb)
    #E0_T_InP = E_T(E0_InP,dt0_InP,b0_InP)
	E0Delta0_T_InAs = E_T(E0Delta0_InAs,dtE0Delta0_InAs,bE0Delta0_InAs)
	E0Delta0_T_InSb = E_T(E0Delta0_InSb,dtE0Delta0_InSb,bE0Delta0_InSb)
	E1_T_InAs = E_T(E1_InAs,dt1_InAs,b1_InAs)
	E1_T_InSb = E_T(E1_InSb,dt1_InSb,b1_InSb)
	E1Delta1_T_InAs = E_T(E1Delta1_InAs,dtE1Delta1_InAs,bE1Delta1_InAs)
	E1Delta1_T_InSb =E_T(E1Delta1_InSb,dtE1Delta1_InSb,bE1Delta1_InSb)        
	E2_T_InAs = E_T(E2_InAs,dt2_InAs,b2_InAs)
	E2_T_InSb = E_T(E2_InSb,dt2_InSb,b2_InSb)
	Eg_T_InAs = E_T(Eg_InAs,dtg_InAs,bg_InAs)
	Eg_T_InSb = E_T(Eg_InSb,dtg_InSb,bg_InSb)
	
	A = func_Q2(x,y,A_InAs,A_InSb,A_InP)
    C= func_Q2(x,y,C_InAs,C_InSb,C_InP)
    D = func_Q2(x,y,D_InAs,D_InSb,D_InP)
    G_E0 = func_Q2(x,y,G_E0_InAs,G_E0_InSb,G_E0_InP) #No value found for Gamma_E1
    G_E1 = func_Q2(x,y,G_E1_InAs,G_E1_InSb,G_E1_InP)
    G_Eg = func_Q2(x,y,G_Eg_InAs,G_Eg_InSb,G_Eg_InP)
    G_E2 = func_Q2(x,y,G_E2_InAs,G_E2_InSb,G_E2_InP)
    E0Delta0 = func_Q2(x,y,E0Delta0_T_InAs, E0Delta0_T_InSb,E0Delta0_InP) #no E0Delta0_T_InP
    E1 = func_Q2(x,y,E1_T_InAs,E1_T_InSb,E1_InP) #no E1_T_InP
    E1Delta1 = func_Q2(x,y,E1Delta1_T_InAs,E1Delta1_T_InSb,E1Delta1_InP)
    Delta1= func_Q2(x,y,Delta1_InAs,Delta1_InSb,Delta1_InP) 
    E2= func_Q2(x,y,E2_T_InAs,E2_T_InSb,E2_InP)
    Eg= func_Q2(x,y,Eg_T_InAs,Eg_T_InSb,Eg_InP)
    alpha_lc = func_Q2(x,y,alpha_lc_InAs,alpha_lc_InSb,alpha_lc_InP)
    =#
    E0 =  0.512 + 0.03*x - 0.183*x^2  #only good for InAsSbP lattice matched to InAs #func_Q2(x,y,E0_T_InAs,E0_T_InSb,E0_T_InP)    
    mstar_ntype = func_Q2(x,y,0.024,0.013,0.07927)*m_0 #effective masses from Adachi, "Properties of Semiconductor alloys: ...," 2009
    mstar_ptype_lh = func_Q2(x,y,0.026,0.014,0.11)*m_0
	mstar_ptype_hh = func_Q2(x,y,0.36,0.38,0.69)*m_0
    mstar_ptype = m_star_ptype(mstar_ptype_hh,mstar_ptype_lh)
    gamma_ntype = Gamma_ntype_InAsSbP(x,y,N0,T,mstar_ntype)
    gamma_ptype = Gamma_ptype_InAsSbP(x,y,N0,T,mstar_ptype)
    epsinf = func_Q2(x,y,eps_inf_InAs,eps_inf_InSb,eps_inf_InP)
    P = func_Q2(x,y,P_InAs,P_InSb,P_InP)
    Nc = 2.0*(3.0*E0*kb*T/(8.0*pi*P^2.0))^(3.0/2.0) #func_Q2(x,y,Nc_InAs,Nc_InSb,Nc_InP)
    #calculating the fermi energy
	res = optimize(t -> fermi(t,E0/(kb*T),pi^(1/2)/2*N0/(Nc*10.0^6)),[E0],Newton()) #Nc from cm-3 to m-3
	F = Optim.minimizer(res)[1]*kb*T + E0 #fermi energy relative to the valence band
    return InAsSbPDsc(x,y,N0,T,E0,gamma_ntype ,gamma_ptype ,mstar_ntype ,mstar_ptype,mstar_ptype_lh,mstar_ptype_hh,epsinf,P,F)
end
precompile(InAsSbP_struct,(Float64,Float64,Float64,Float64,))

@inline function eVtoOmega(energy) #energy in eV
    return(energy*evJ/hbar)
end
precompile(eVtoOmega, (Float64,))



@inline function  eps_InAsSbP_xy_ntype(E_photon,InAsSbPstruct)

    #adding IB or no IB. IB - epsinf since I don't want to double count
    #if(E_photon>InAsSbPstruct.E0_InAsSbP)
    return epsIB(eVtoOmega(E_photon),InAsSbPstruct.N0,InAsSbPstruct.T,InAsSbPstruct.E0_InAsSbP,InAsSbPstruct.epsinf,InAsSbPstruct.P,InAsSbPstruct.mstar_ptype_hh,InAsSbPstruct.mstar_ptype_lh,InAsSbPstruct.F,0.0) - InAsSbPstruct.epsinf + epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(E_photon),InAsSbPstruct.mstar_ntype,InAsSbPstruct.gamma_ntype,InAsSbPstruct.N0,InAsSbPstruct.epsinf)
    #else
    #    return epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(E_photon),InAsSbPstruct.mstar_ntype,InAsSbPstruct.gamma_ntype,InAsSbPstruct.N0,InAsSbPstruct.epsinf)
    #end
end
precompile(eps_InAsSbP_xy_ntype,(Float64,InAsSbPDsc))


@inline function  eps_InAsSbP_xy_ptype(E_photon,InAsSbPstruct)
    
    #adding the IB and FCL to eps
    #if(E_photon>InAsSbPstruct.E0_InAsSbP)
    return epsIB(eVtoOmega(E_photon),InAsSbPstruct.N0,InAsSbPstruct.T,InAsSbPstruct.E0_InAsSbP,InAsSbPstruct.epsinf,InAsSbPstruct.P,InAsSbPstruct.mstar_ptype_hh,InAsSbPstruct.mstar_ptype_lh,InAsSbPstruct.F,1.0) - InAsSbPstruct.epsinf + epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(E_photon),InAsSbPstruct.mstar_ptype,InAsSbPstruct.gamma_ptype,InAsSbPstruct.N0,InAsSbPstruct.epsinf)
    #else
    #    return epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(E_photon),InAsSbPstruct.mstar_ptype,InAsSbPstruct.gamma_ptype,InAsSbPstruct.N0,InAsSbPstruct.epsinf)
    #end
end
precompile(eps_InAsSbP_xy_ptype,(Float64,InAsSbPDsc))



@inline function eps_InAsSbP_imag_xy_ptype(enr,InAsSbPstruct)
    #In the code InAs_xSb_yP_1-x-y,  from III-V Ternary and Quaternary Compounds Table 30.14 uses InP_xAs_ySb_1-x-y, and E0 = 0.512+0.030*y-0.183*y^2
    #E_bandgap = 0.512 + 0.030*InAsSbPstruct.x-0.183*InAsSbPstruct.x^2   
    if(enr>InAsSbPstruct.E0_InAsSbP)
        return imag(epsIB(eVtoOmega(enr),InAsSbPstruct.N0,InAsSbPstruct.T,InAsSbPstruct.E0_InAsSbP,InAsSbPstruct.epsinf,InAsSbPstruct.P,InAsSbPstruct.mstar_ptype_hh,InAsSbPstruct.mstar_ptype_lh,InAsSbPstruct.F,1.0))
    else
        return 0.0
    end
end
precompile(eps_InAsSbP_imag_xy_ptype,(Float64,InAsSbPDsc))




@inline function eps_InAsSbP_imag_xy_ntype(enr,InAsSbPstruct)
    #In the code InAs_xSb_yP_1-x-y,  from III-V Ternary and Quaternary Compounds Table 30.14 uses InP_xAs_ySb_1-x-y, and E0 = 0.512+0.030*y-0.183*y^2
    #E_bandgap = 0.512 + 0.030*InAsSbPstruct.x-0.183*InAsSbPstruct.x^2   
    if(enr>InAsSbPstruct.E0_InAsSbP)
        return imag(epsIB(eVtoOmega(enr),InAsSbPstruct.N0,InAsSbPstruct.T,InAsSbPstruct.E0_InAsSbP,InAsSbPstruct.epsinf,InAsSbPstruct.P,InAsSbPstruct.mstar_ptype_hh,InAsSbPstruct.mstar_ptype_lh,InAsSbPstruct.F,0.0))
    else
        return 0.0
    end
end
precompile(eps_InAsSbP_imag_xy_ntype,(Float64,InAsSbPDsc))



@inline function epsFCL_InAsSbP(x,y,omega,mstar,gamma,N_base,eps_inf)
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
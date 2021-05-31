__precompile__()
module eps_InAsSbP_v2
export func_Q,func_Q2,E_T,func_f,InAsSbP_struct, eps_InAsSbP_xy_ntype,eps_InAsSbP_xy_ptype,eps_InAsSbP_imag_xy,Gamma_ntype_InAsSbP,Gamma_ptype_InAsSbP,eVtoOmega,epsFCL_InAsSbP
#Used this version in the end, accurate imaginary part from the range of 0.1 to 1.5, above 1.5, imaginary part becomes negative"""
#version 2: Includes free carrier and lattice absorption. Also used InP parameters from Gavin's fit.
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

#list of constants: list of constants for InAsSbP, AD=InAs, BD=InSb, CD=InP ABC=AsSbP  #Cuevas Paper Table1
#all parameters from fits to binary materials

epsinf_InAs = 0.489
epsinf_InSb = 0.0
epsinf_InP = 1.334

A_InAs = 62.787   #A (eV^1.5) for eps_1, Cuevas
A_InSb = 4.121   #A (eV^1.5) for eps_1
A_InP =  4.738    # Fit to Adachi data

C_InAs = 1.501  #C(Dimension-less) for eps_3 Cuevas
C_InSb = 3.994
C_InP=  1.724    #Fit to Adachi data

D_InAs = 543.233 #D(Dimension-less) for eps_4 Cuevas
D_InSb = 83.671 
D_InP=  102.135   #Fit to Adachi data

G_E0_InAs = 96.385 #in eV Cuevas
G_E0_InSb = 706.646
G_E0_InP= 0.017 #Fit to Adachi data

G_E1_InAs = 0.074 #
G_E1_InSb = 0.193 
G_E1_InP  = 0.005 # 

G_E2_InAs =  0.442 #gamma in eV 
G_E2_InSb =  1.218  
G_E2_InP =  0.509 #

G_Eg_InAs =  4.388 #(eV) 
G_Eg_InSb = 1.282 
G_Eg_InP  = 0.113

B1_InAs = 4.360
B1_InSb = 8.747
B1_InP = 4.28

B2_InAs = 0.0
B2_InSb = 0.009
B2_InP = 0.0

#bandgap as input parameters to the fit, with ternary bowing constants from Adachi, "Properties of semiconductor alloys"

#Cuevas Table2, The energies are in eV, NO bowing constants for E0??Only temperature depence?
E0_InAs = 0.42  #at 300K, E0_T = 1.5,at 0K 1.4
dt0_InAs = 0.28*0.001
b0_InAs = 93.0
E0_InSb = 0.24
dt0_InSb = 0.32*0.001
b0_InSb = 170.0  #Need fitting Adachi Table1
E0_InP = 1.35  
E0_InAsP = 0.145

E0Delta0_InAs = 0.79#-0.42 
dtE0Delta0_InAs = 0.34*0.001
bE0Delta0_InAs  = 248.0
E0Delta0_InSb = 1.20#-0.24
dtE0Delta0_InSb = 0.32*0.001
bE0Delta0_InSb = 170.0 #E0Delta0_GaAs = 1.85#-0.152 dtE0Delta0_GaAs = 0.35*0.001 bE0Delta0_GaAs = 225.0 E0Delta0_GaSb = 1.57#-0.81 dtE0Delta0_GaSb = 0.42*0.001
E0Delta0_InP = 1.45 # Adachi Table B15-1
E0Delta0_InAsP = 0.145


E1_InAs = 2.61 #Cuevas Table2 
dt1_InAs = 0.5*0.001
b1_InAs = 0.0
E1_InSb = 2.0 
dt1_InSb = 0.68*0.001
b1_InSb = 132.0
E1_InP = 3.17 # Adachi Table B15-1
E1_InAsP = 0.25
  
E1Delta1_InAs =2.88#-2.61 #Cuevas Table 2
dtE1Delta1_InAs = 0.5*0.001
bE1Delta1_InAs = 0.0
E1Delta1_InSb = 2.49#-2.00
dtE1Delta1_InSb = 0.65*100
bE1Delta1_InSb = 170.0 
E1Delta1_InP = 3.29 #Adachi  Table B15-1
E1Delta1_InAsP = 0.25
  
Delta1_InAs = 2.73-2.46 #Cuevas
Delta1_InSb = 2.858-2.366 #Delta1_InP = 0.15 #Adachi  Table1
Delta1_InP =  3.29-3.17 
Delta1_InAsP = 0.25

E2_InAs = 4.74 #from Zhang Thermo Conversion
dt2_InAs = 0.56*0.001
b2_InAs = 0.0
E2_InSb = 4.24
dt2_InSb = 0.54*0.001
b2_InSb = 0.0
E2_InP = 4.7 #Adachi  Table B15-1
E2_InAsP = 0.03

Eg_InAs = 1.13 #Cuevas Table 2
dtg_InAs = 0.28*0.001
bg_InAs = 93.0
Eg_InSb = 0.93
dtg_InSb = 0.32*0.001
bg_InSb = 170.0
Eg_InP  = 2.05 #Adachi Table B15-1
Eg_InAsP = 0.145

#Cuevas paper Table IV in angstrong
alpha_lc_InAs = 6.06
alpha_lc_InSb = 6.48
alpha_lc_InP =  5.8697 
#alpha_lc_GaAs = 5.65 alpha_lc_GaSb = 6.10 alpha_lc_GaIn = 0.055  #Cuevas no value #zhang 0.055 alpha_lc_AsSb = -0.023 #Cuevas no value #zhang -0.023



#Q Parameter in Notes, Using Ternary bowing params, AB=InAs, AC=InSb, AD=InP BC=AsSb, BD=AsP, CD=SbP
@inline function  func_Q(x,y,B_AB,B_AC,B_AD,C_BC,C_BD,C_CD)
    return x*B_AB + y*B_AC + (1-x-y)*B_AD  + C_BC*x*y + C_BD*x*(1-x-y) + C_CD*y*(1-x-y)
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

@inline function Gamma_ntype_InAsSbP(x,y,N_Base,T)
    
    #InAs
    me_InAs = 0.023*m_0
    umin_InAs = 1000.0
    umax_InAs = 34000.0
    Nref_InAs = 1.1*10^18*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InAs = 0.32
    t1_InAs = 1.57
    t2_InAs = 3.0

    #InSbP (same as Sentaurus)
    me_InSbP = 0.063*m_0
    umin_InSbP = 559.6
    umax_InSbP = 6756.9
    Nref_InSbP = 2.287*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InSbP = 0.5413
    t1_InSbP = 2.0
    t2_InSbP = 3.1105

    #interpolated values
    m_star = x*me_InAs + (1-x)*me_InSbP
    umin = (x/umin_InAs + (1-x)/umin_InSbP)
    umax = (x/umax_InAs + (1-x)/umax_InSbP)
    Nref = x*Nref_InAs + (1-x)*Nref_InSbP
    phi = x*phi_InAs + (1-x)*phi_InSbP
    t1 = x*t1_InAs + (1-x)*t1_InSbP
    t2 = x*t2_InAs + (1-x)*t2_InSbP

    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    return e/(m_star*mu_franc*10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ntype_InAsSbP, (Float64,Float64,Float64,Float64))

@inline function Gamma_ptype_InAsSbP(x,y,N_Base,T)
    
    #InAs
    me_InAs = 0.41*m_0
    umin_InAs = 1000.0
    umax_InAs = 34000.0
    Nref_InAs = 1.1*10^18*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InAs = 0.32
    t1_InAs = 1.57
    t2_InAs = 3.0

    #InSbP (same as Sentaurus)
    me_InSbP = 0.54*m_0
    umin_InSbP = 13.87
    umax_InSbP = 208.8
    Nref_InSbP = 5.2203*10^17*10^6 #in m^(-3) since N_base in m^(-3)
    phi_InSbP = 0.6138
    t1_InSbP = 1.907
    t2_InSbP = 2.907

    #interpolated values
    m_star = x*me_InAs + (1-x)*me_InSbP
    umin = (x/umin_InAs + (1-x)/umin_InSbP)
    umax = (x/umax_InAs + (1-x)/umax_InSbP)
    Nref = x*Nref_InAs + (1-x)*Nref_InSbP
    phi = x*phi_InAs + (1-x)*phi_InSbP
    t1 = x*t1_InAs + (1-x)*t1_InSbP
    t2 = x*t2_InAs + (1-x)*t2_InSbP

    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    return e/(m_star*mu_franc*10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ptype_InAsSbP, (Float64,Float64,Float64,Float64))

struct InAsSbPDsc
    x::Float64
    y::Float64
    N0::Float64
    T::Float64
    A_InAsSbP::Float64
    C_InAsSbP::Float64
    D_InAsSbP::Float64
    B1_InAsSbP::Float64
    B2_InAsSbP::Float64
    G_E0_InAsSbP::Float64 #No value found for Gamma_E1
    G_E1_InAsSbP::Float64
    G_Eg_InAsSbP::Float64
    G_E2_InAsSbP::Float64
    E0_InAsSbP::Float64
    E0Delta0_InAsSbP::Float64
    E1_InAsSbP::Float64
    E1Delta1_InAsSbP::Float64
    Delta1_InAsSbP::Float64
    E2_InAsSbP::Float64
    Eg_InAsSbP::Float64
    alpha_lc_InAsSbP::Float64
    gamma_ntype::Float64
    gamma_ptype::Float64
    mstar_ntype::Float64
    mstar_ptype::Float64
    epsinf::Float64
end

@inline function InAsSbP_struct(x,y,N0,T)
     #assume equal proportion of As,Sb and P
    #list of constants for InAsSbP, AD=InAs, BD=InSb, CD=InP ABC=AsSbP  #Cuevas Paper Table1
    #Bowing Constants for GaIn, AsSb from ZhangThermoConversion
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
    B1 = func_Q2(x,y,B1_InAs,B1_InSb,B1_InP)
    B2 = func_Q2(x,y,B2_InAs,B2_InSb,B2_InP)
    G_E0 = func_Q2(x,y,G_E0_InAs,G_E0_InSb,G_E0_InP) 
    G_E1 = func_Q2(x,y,G_E1_InAs,G_E1_InSb,G_E1_InP)
    G_Eg = func_Q2(x,y,G_Eg_InAs,G_Eg_InSb,G_Eg_InP)
    G_E2 = func_Q2(x,y,G_E2_InAs,G_E2_InSb,G_E2_InP)
    E0 =  0.512 + 0.03*x - 0.183*x^2  #only good for InAs substrate #func_Q2(x,y,E0_T_InAs,E0_T_InSb,E0_T_InP)    
    E0Delta0 = func_Q2(x,y,E0Delta0_T_InAs, E0Delta0_T_InSb,E0Delta0_InP) 
    E1 = func_Q2(x,y,E1_T_InAs,E1_T_InSb,E1_InP) 
    E1Delta1 = func_Q2(x,y,E1Delta1_T_InAs,E1Delta1_T_InSb,E1Delta1_InP)
    Delta1= func_Q2(x,y,Delta1_InAs,Delta1_InSb,Delta1_InP) 
    E2= func_Q2(x,y,E2_T_InAs,E2_T_InSb,E2_InP)
    Eg= func_Q2(x,y,Eg_T_InAs,Eg_T_InSb,Eg_InP)
    alpha_lc = func_Q2(x,y,alpha_lc_InAs,alpha_lc_InSb,alpha_lc_InP)
    gamma_ntype = Gamma_ptype_InAsSbP(x,y,N0,T)
    gamma_ptype = Gamma_ptype_InAsSbP(x,y,N0,T)
    mstar_ntype = x*0.023*m_0 + (1-x)*0.063*m_0
    mstar_ptype = x*0.41*m_0 + (1-x)*0.54*m_0
    epsinf = func_Q2(x,y,epsinf_InAs,epsinf_InSb,epsinf_InP)
    return InAsSbPDsc(x,y,N0,T,A,C,D,B1,B2,G_E0,G_E1,G_Eg,G_E2,E0,E0Delta0,E1,E1Delta1, Delta1, E2,Eg ,alpha_lc ,gamma_ntype ,gamma_ptype ,mstar_ntype ,mstar_ptype,epsinf)
end
precompile(InAsSbP_struct,(Float64,Float64,Float64,Float64))

@inline function eVtoOmega(energy) #energy in eV
    return(energy*evJ/hbar)
end
precompile(eVtoOmega, (Float64,))

@inline function  eps_InAsSbP_xy(E_photon,InAsSbPstruct)
    #GaInAsSb in the energy range 0.5-6eV omega  (hbar*2*math.pi)*c/(0.5*eV) = 2.5um, frequency 

    eps_1 =InAsSbPstruct.A_InAsSbP*InAsSbPstruct.E0_InAsSbP^(-1.5)*(func_f((E_photon + im*InAsSbPstruct.G_E0_InAsSbP)/InAsSbPstruct.E0_InAsSbP) + 0.5*(InAsSbPstruct.E0_InAsSbP/(InAsSbPstruct.E0Delta0_InAsSbP))^1.5*func_f((E_photon + im*InAsSbPstruct.G_E0_InAsSbP)/InAsSbPstruct.E0Delta0_InAsSbP))
    #B1=44*(InAsSbPstruct.E1_InAsSbP+InAsSbPstruct.Delta1_InAsSbP/3)/(InAsSbPstruct.alpha_lc_InAsSbP*InAsSbPstruct.E1_InAsSbP^2)
    #B2=44*(InAsSbPstruct.E1_InAsSbP+2*InAsSbPstruct.Delta1_InAsSbP/3)/(InAsSbPstruct.alpha_lc_InAsSbP*(InAsSbPstruct.E1_InAsSbP+InAsSbPstruct.Delta1_InAsSbP)^2)
    eps_2term1 = -InAsSbPstruct.B1_InAsSbP*((E_photon + im*InAsSbPstruct.G_E1_InAsSbP)/InAsSbPstruct.E1_InAsSbP)^(-2)*log(1-((E_photon + im*InAsSbPstruct.G_E1_InAsSbP)/InAsSbPstruct.E1_InAsSbP)^2)
    eps_2term2 = -InAsSbPstruct.B2_InAsSbP*((E_photon + im*InAsSbPstruct.G_E1_InAsSbP)/(InAsSbPstruct.E1Delta1_InAsSbP))^(-2)*log(1-((E_photon + im*InAsSbPstruct.G_E1_InAsSbP)/(InAsSbPstruct.E1Delta1_InAsSbP))^2)
    eps_2 = eps_2term1 + eps_2term2 
    eps_3 = InAsSbPstruct.C_InAsSbP*InAsSbPstruct.E2_InAsSbP^2/(InAsSbPstruct.E2_InAsSbP^2-E_photon^2 + im*(-E_photon*InAsSbPstruct.G_E2_InAsSbP))
    E_complex = (E_photon + im*InAsSbPstruct.G_Eg_InAsSbP)
    eps_4_1st = -log(InAsSbPstruct.E1_InAsSbP/InAsSbPstruct.Eg_InAsSbP)*InAsSbPstruct.Eg_InAsSbP^2/E_complex^2
    eps_4_2nd = (1+InAsSbPstruct.Eg_InAsSbP/E_complex)^2*log((E_complex+InAsSbPstruct.E1_InAsSbP)/(E_complex+InAsSbPstruct.Eg_InAsSbP))/2
    eps_4_3rd = (1-InAsSbPstruct.Eg_InAsSbP/E_complex)^2*log((E_complex-InAsSbPstruct.E1_InAsSbP)/(E_complex-InAsSbPstruct.Eg_InAsSbP))/2
    eps_4 = (2*InAsSbPstruct.D_InAsSbP/pi)*(eps_4_1st+eps_4_2nd+eps_4_3rd)
    
    eps = InAsSbPstruct.epsinf + eps_1 + eps_2 + eps_3 + eps_4 #+ epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(E_photon),InAsSbPstruct.mstar_ntype,InAsSbPstruct.gamma_ntype,InAsSbPstruct.N0)
    return eps
end
precompile(eps_InAsSbP_xy,(Float64,InAsSbPDsc))


@inline function  eps_InAsSbP_xy_ptype(enr,InAsSbPstruct)
    if(enr>InAsSbPstruct.E0_InAsSbP)
        return eps_InAsSbP_xy(enr,InAsSbPstruct) + epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(enr),InAsSbPstruct.mstar_ptype,InAsSbPstruct.gamma_ptype,InAsSbPstruct.N0)
    else
        return real(eps_InAsSbP_xy(enr,InAsSbPstruct)) + epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(enr),InAsSbPstruct.mstar_ptype,InAsSbPstruct.gamma_ptype,InAsSbPstruct.N0)
    end
end
precompile(eps_InAsSbP_xy_ptype,(Float64,InAsSbPDsc))

@inline function  eps_InAsSbP_xy_ntype(enr,InAsSbPstruct)
    if(enr>InAsSbPstruct.E0_InAsSbP)
        return eps_InAsSbP_xy(enr,InAsSbPstruct) + epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(enr),InAsSbPstruct.mstar_ntype,InAsSbPstruct.gamma_ntype,InAsSbPstruct.N0)
    else
        return real(eps_InAsSbP_xy(enr,InAsSbPstruct)) + epsFCL_InAsSbP(InAsSbPstruct.x,InAsSbPstruct.y,eVtoOmega(enr),InAsSbPstruct.mstar_ntype,InAsSbPstruct.gamma_ntype,InAsSbPstruct.N0)
    end
end
precompile(eps_InAsSbP_xy_ntype,(Float64,InAsSbPDsc))

@inline function eps_InAsSbP_imag_xy(enr,InAsSbPstruct)
    #In the code InAs_xSb_yP_1-x-y,  from III-V Ternary and Quaternary Compounds Table 30.14 uses InP_xAs_ySb_1-x-y, and E0 = 0.512+0.030*y-0.183*y^2
    #E_bandgap = 0.512 + 0.030*InAsSbPstruct.x-0.183*InAsSbPstruct.x^2   
    if(enr>InAsSbPstruct.E0_InAsSbP)
        return imag(eps_InAsSbP_xy_ntype(enr,InAsSbPstruct))
    else
        return 0
    end
end
precompile(eps_InAsSbP_imag_xy,(Float64,InAsSbPDsc))




@inline function epsFCL_InAsSbP(x,y,omega,mstar,gamma,N_base)
    
    #eps,gamma parameters from Sadao Adachi: Optical Properties of Crystalline and Amorphous Semiconductors (1999)
    #To and LO parameters: Adachi - Properties of semiconductor alloys
    #InAs
    eps_inf_InAs = 12.3
    g_InAs = 9.24*10.0^11  #s^(-1), gamma
    o_TO_InAs = 4.14*10.0^(13) #s^(-1)
    o_LO_InAs = 4.55*10.0^(13)


    #InSb
    eps_inf_InSb = 15.68
    g_InSb = 5.65*10.0^11  #s^(-1), gamma
    o_TO_InSb = 3.39*10.0^(13) #s^(-1)
    o_LO_InSb =3.59*10.0^(13)


    #InP
    eps_inf_InP = 9.66
    g_InP = 3.58*10.0^11  #s^(-1), gamma
    o_TO_InP = 5.74*10.0^(13) #s^(-1)
    o_LO_InP = 6.53*10.0^(13)
    
    eps_inf = x*eps_inf_InAs + y*eps_inf_InSb + (1-x-y)*eps_inf_InP 
    o_p_square=N_base*e^2/(eps_0*eps_inf*mstar) #s^(-2)

    return - o_p_square/(omega*(omega+ im*gamma)) + x*eps_inf_InAs*(o_LO_InAs^2-o_TO_InAs^2)/(o_TO_InAs^2-omega^2-im*omega*g_InAs)+y*eps_inf_InSb*(o_LO_InSb^2-o_TO_InSb^2)/(o_TO_InSb^2-omega^2-im*omega*g_InSb)+(1-x-y)*eps_inf_InP*(o_LO_InP^2-o_TO_InP^2)/(o_TO_InP^2-omega^2-im*omega*g_InP)
end
precompile(epsFCL_InAsSbP,(Float64,Float64,Float64,Float64,Float64,Float64))
end
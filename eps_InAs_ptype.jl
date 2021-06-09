__precompile__()
module eps_InAs_ptype
using QuadGK
export epsFCL_ptype,Gamma_ptype,m_star_ptype,E0_T_InAs_ptype
#from integration import Integrate
#from scipy.integrate import quad
 #but I can't use QuadGk for larger values of eV, can only functionine function Integrate, tested E = 0.0037eV, get same eps of 10.240628772514803 + 0.03160036047883907im

const e = 1.602*^(10,-19)
const eV = 1.60218*^(10,-19)
const eps_0 = 8.8542*^(10,-12)
const c = 3*^(10,8)
const hbar = 1.05457*^(10,-34) #m2kg/s
const hbar_eV = 6.5821*^(10,-16) #eV s 
const m_e =9.109*^(10,-31)#kg
const m_0 = 0.9109*^(10,-30)
const m_0_eV = 0.51099906*^(10,6) #eV

const E0_InAs = 0.42
const dt0_InAs = 0.28*0.001
const b0_InAs = 93
const A_InAs = 8.34*10^(-9) #for p-type  InAs, From Levinshtein Handbook on semiconductor parameters.
const B_InAs = 2.91*10^(-7)
const C_InAs = 4.53*10^(-12) #N in m^(-3)

const P_square = 11.9   #eV
const E_g_InAs= 0.416 #eV 
const m_h = 0.41      #effective mass of holes 

const umin =20.0
const umax = 530.0
const Nref = 1.1*10^17*10^6 #in m^(-3) since N_base in m^(-3)
const phi = 0.46
const t1 = 2.3
const t2 = 3.0


@inline function m_star_ptype(m_hh,m_lh) # change m_e according to MB shift
    #term = 8.0*P_square/eV*hbar^2*(3.0*pi^2*n)^(2/3)/(3.0*m_e*E_g_InAs^2) 
    return    (m_hh^(3/2)+m_lh^(3/2))^(2/3) #(1.0+(4.0*P_square/(3.0*E_g_InAs))*(1.0+term)^(-1/2))^(-1)*m_0
end 
precompile(m_star_ptype, (Float64,Float64))

@inline function Gamma_ptype(N_Base,T,m_star_ptype)
    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    return e/(m_star_ptype*mu_franc*10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ptype, (Float64,Float64,Float64))

@inline function E0_T_InAs_ptype(N0,T) 
    #delta_Eg_InAs = A_InAs*(N0/10^6)^(1/3) + B_InAs*(N0/10^6)^(1/4) + C_InAs*(N0/10^6)^(1/2)

    #We don't include the effect of Moss Burstein shift for p-type since high VB density of states
    #dont include bandgap narrowing, since the states would not contribute to optical transistions
    E_T = E0_InAs - dt0_InAs*T^2/(T+b0_InAs)
    return E_T #- delta_Eg_InAs      #changed from + to - since it is bandgap narrowing 
end
precompile(E0_T_InAs_ptype, (Float64,Float64))

const eps_inf_InAs = 12.3
const g = 2.89*10.0^11  #s^(-1), gamma
const o_TO = 4.12*10.0^(13) #s^(-1)
const o_LO = 4.58*10.0^(13)
@inline function epsFCL_ptype(omega,mstar,gamma,N_base)
    o_p_square=N_base*e^2/(12.3*eps_0*mstar) #s^(-2)
    return eps_inf_InAs*(1.0+(o_LO^2-o_TO^2)/(o_TO^2-omega^2-im*omega*g)-o_p_square/(omega*(omega+ im*gamma)))
end
precompile(epsFCL_ptype,(Float64,Float64,Float64,Float64))
end
#@inline function epsIB(omega,N0,T,E0_T_InAs_value) #use the function from ntype
#    if(hbar_eV*omega<=E0_T_InAs_value)
#        alpha_lb_2=0.0
#    else
#        alpha_lb_2= (3*10^8*(hbar_eV*omega-E0_T_InAs_value))^(1/2)*100
#    end
#    m_IB_pp = alpha_lb_2*c/(2*omega)

#    integratedfn(omega_p) = alpha_lb_2(omega_p,N0,T)/(omega^2-omega_p^2)
#    if(hbar_eV*omega<E0_T_InAs_value)
#        m_IB_p_value =1.0
#    else
     #   result  =1.0 + (c/pi)*quadgk(integratedfn,100.0,10.0^25,rtol=1e-8)[1] # result  =(1.0 + (c/pi)*Integrate(integratedfn, 10000,100.0,10.0^25))
#        m_IB_p_value = (1.0 + 0.00026172525)
#    end
#    return (m_IB_p_value^2-m_IB_pp^2) + im*2*m_IB_p_value*m_IB_pp
#end
#precompile(epsIB,(Float64,Float64,Float64,Float64))

#function eps_InAs_ntype(E_photon::Float64,N_base::Float64)::ComplexF64 #N_base in m^(-3)
    #E_photon= eVtoOmega(omega)
#    omega = eVtoOmega(E_photon)
#    T=300.0 #o_angular = 5.38*10^14   #hbar*5.38*10^14/eV = 0.35eV
#    if(E_photon>0.35)
#        return epsIB(omega,N_base,T) +  epsFCL(omega,N_base,T)
#    else
#        return epsFCL(omega,N_base,T)
#    end
#end

#omegalist= 2.0*10.0^14:5.75*10.0^13:2.5*10.0^15#
#eVlist=  0.2:0.02:1.5
#epslist = []
#for ind = 1:length(eVlist)
#    eps = eps_InAs_ntype(eVlist[ind],10.0^23)
#    println("eps_InAs",eps)
#    push!(epslist, imag(eps))
#end
#plot(eVlist,epslist)
#savefig("test_InAs_real.png")
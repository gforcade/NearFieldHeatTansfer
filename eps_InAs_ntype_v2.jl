__precompile__()
module eps_InAs_ntype_v2
using QuadGK,Cubature
export OmegatoeV, m_star, Gamma_ntype, E0_T_InAs,epsFCL,epsIB,epsIBEV,fermi
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
const kb = 8.617333*^(10,-5) #eV K-1

const E0_InAs = 0.417
const dt0_InAs = 0.276*0.001
const b0_InAs = 93
const A_InAs = 4*10^(-8)  #for n-type  InAs, N in cm^(-3)
const B_InAs = 0.0 #1.97*10^(-7)
const C_InAs = 0.0# 57.9*10^(-12)

@inline function OmegatoeV(omega)
    return omega*hbar_eV
end
precompile(OmegatoeV,(Float64,))

@inline function m_star(n) # change m_e according to MB shift
    P_square = 11.9   #eV
    E_g= 0.416 #eV 
    term = 8.0*P_square/eV*hbar^2*(3.0*pi^2*n)^(2/3)/(3.0*m_e*E_g^2) 
    return (1.0+(4.0*P_square/(3.0*E_g))*(1.0+term)^(-1/2))^(-1)*m_0
end 
precompile(m_star, (Float64,))

@inline function Gamma_ntype(N_Base,T)
    umin = 1000.0
    umax = 34000.0
    Nref = 1.1*10^18*10^6 #in m^(-3) since N_base in m^(-3)
    phi = 0.32
    t1 = 1.57
    t2 = 3.0
    mu_franc =  umin + (umax*(300.0/T)^t1-umin)/(1+(N_Base/(Nref*(T/300.0)^t2))^phi)#from Francoeur 
    return e/(m_star(N_Base)*mu_franc*10000.0) #s^(-1) change mu from cm^(-2) to m^(-2) 
end
precompile(Gamma_ntype, (Float64,Float64))

@inline function E0_T_InAs(N0,T) 
    #bandgap narrowing #Gavin: Got rid of bandgap narrowing and moss burstein shift to have better fit with experimental data
    #delta_Eg_InAs = A_InAs*(N0/10^6)^(1/3) + B_InAs*(N0/10^6)^(1/4) + C_InAs*(N0/10^6)^(1/2)
    #mstar = m_star(N0)
    #Need to include the effect of Moss Burstein shift on E0_InAs
    E_T = E0_InAs -dt0_InAs*T^2/(T+b0_InAs)  #   + (hbar^2/(2*mstar))*(3*N0/(8*pi))^(2/3) #changed pi^2 to 1/8pi (was an error)
    return E_T #- delta_Eg_InAs     #changed from + to - since it is bandgap narrowing 
end
precompile(E0_T_InAs, (Float64,Float64))



function fermi(t,phi,N)
    #calculates the n-type fermi energy for nonparabolic conduction band, with alpha = 1/E_0
    #integration function
    func(x) = x^(1/2)*(1.0 + x/phi)^(1/2)*(1.0 + 2.0*x/phi)/(1.0 + exp(x - t[1]))
    integ = quadgk(func,0.0,10.0^2)[1]
    return abs(integ - N)       
end
precompile(fermi,(Float64,Float64,Float64))



@inline function k_omega(omega,E0_T_InAs_value,P,mh)
    #equation in Anderson 1980 paper (Eq. 18)
    #breaking down equation  #divide P by 100 to convert from eVcm to eVm #divide P term by e to convert to eV3 s2 kg-1 
    lump = (1.0 + m_0/mh)*(2*hbar_eV*omega/E0_T_InAs_value - 1.0)
    k1 = (4.0/3.0*(P/100.0)^2/e + hbar_eV^2*E0_T_InAs_value/m_0*lump)/ (hbar_eV^4/m_0^2 * (1.0 + m_0/mh)^2)
    k2 = 1.0 - (1.0 - (4.0*hbar_eV^4/(m_0)^2*(1 + m_0/mh)^2*hbar_eV*omega*(hbar_eV*omega - E0_T_InAs_value))/(4.0/3.0*(P/100.0)^2/e + hbar_eV^2*E0_T_InAs_value/m_0*lump)^2 )^(1/2)
    return k1*k2  #outputs k^2 [kg eV-1 s-2] 
end
precompile(k_omega,(Float64,Float64,Float64,Float64))




@inline function BM(omega,T,mh,k_om,F)
    #Burstein-Moss effect #Anderson 1980 Eq. 17
    top = -expm1(-hbar_eV*omega/(kb*T))
    bot = 1 + exp(- (F + hbar_eV^2*k_om/(2*mh))/(kb*T))
    bot2 = 1 + exp(- (hbar_eV*omega - F - hbar_eV^2*k_om/(2*mh))/(kb*T) )
    return top/(bot*bot2)
end
precompile(BM,(Float64,Float64,Float64,Float64,Float64))



@inline function epsFCL(omega,mstar,gamma,N_base)
    #InAs parameters
    eps_inf_InAs = 11.6 
    g = 9.23*10.0^11  #s^(-1), gamma #from Adachi, "Optical Constants of...", 1999
    o_TO = 4.14*10.0^(13) #s^(-1) #0.0271 [eV] #from Adachi, "Properties of semiconductor and their alloys...", 2009
    o_LO = 4.55*10.0^(13) #0.0301 [eV]
    o_p_square=N_base*e^2/(eps_0*eps_inf_InAs*mstar) #s^(-2)
    return eps_inf_InAs*(1.0-o_p_square/(omega*(omega+ im*gamma))+(o_LO^2-o_TO^2)/(o_TO^2-omega^2-im*omega*g))
end
precompile(epsFCL,(Float64,Float64,Float64,Float64))


@inline function alpha_lh(omega,E0_T_InAs_value,eps_inf,P)
    #function for light hole band absorption #Anderson 1980 Eq. 14
    return (1 + 2*(E0_T_InAs_value/(hbar_eV*omega))^2)*((hbar_eV*omega)^2 - E0_T_InAs_value^2)^(1/2)/(137*(6*eps_inf)^(1/2)*4*P)*100  #multiplied by 100 to convert cm to m
end
precompile(alpha_lh,(Float64,Float64,Float64,Float64))



@inline function alpha_hh(omega,E0_T_InAs_value,eps_inf,P,mhh,k_om)
    #function for heavy hole band absorption #Anderson 1980 Eq. 16
    #the approximate version
    return (3/2*hbar_eV*omega*(hbar_eV*omega - E0_T_InAs_value))^(1/2)/(137*eps_inf^(1/2)*P*(1 + 3/(4*m_0*P^2)*e*100^2*hbar_eV^2*E0_T_InAs_value*(1 + m_0/mhh)*(2*hbar_eV*omega/E0_T_InAs_value - 1)))*100 #multiply by 100 to convert cm to m, and by e to convert eV to J
    #the exact version 
    #sqr = (1 + 8.0*(P/100.0)^2*(k_om/e)/(3.0*E0_T_InAs_value^2))^(1/2)
    #top = E0_T_InAs_value*(k_om/e)^(1/2)/(hbar_eV*omega*2.0)*(sqr + 1) #k from kg eV-1 s-2 to m-1
    #bot = 137*eps_inf^(1/2)*(1.0 + 3.0/4.0*hbar_eV^2*E0_T_InAs_value*e*(1 + m_0/mhh)*sqr/(m_0*(P/100.0)^2))
    #return top/bot 
end
precompile(alpha_hh,(Float64,Float64,Float64,Float64,Float64,Float64))



@inline function alpha_lb_2(omega,E0_T_InAs_value,eps_inf,P,mhh,mlh,T,F,n_or_p)
    #total absorption coefficienct for IB
    if(hbar_eV*omega<=E0_T_InAs_value)
        #below bandgap, no absorption
        alpha_lb=0.0
    else
        #calculate k_omega^2 with units kg eV-1 s-2
        k_om = k_omega(omega,E0_T_InAs_value,P,mhh)
        if n_or_p == 0.0    #using moss-burstein effect when n_type
            alpha_lb= alpha_lh(omega,E0_T_InAs_value,eps_inf,P)*BM(omega,T,mlh,k_om,F) + alpha_hh(omega,E0_T_InAs_value,eps_inf,P,mhh,k_om)*BM(omega,T,mhh,k_om,F)   #(3*10^8*(hbar_eV*omega-E0_T_InAs_value))^(1/2)*100 , legacy version
        else #p-type material, thus dont include moss-burstein effect
            alpha_lb= alpha_lh(omega,E0_T_InAs_value,eps_inf,P) + alpha_hh(omega,E0_T_InAs_value,eps_inf,P,mhh,k_om)
        end
    end
    return alpha_lb
end
precompile(alpha_lb_2,(Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64))



function epsIB(omega,N0,T,E0_T_InAs_value,eps_inf,P,mhh,mlh,F,n_or_p)
    #dielectric function from interband transitions
    # - n_or_p states if the material is n type (0.0) or p type (1.0)
    #extinction coefficient
    alpha = alpha_lb_2(omega,E0_T_InAs_value,eps_inf,P,mhh,mlh,T,F,n_or_p)
    m_IB_pp = alpha*c/(2*omega)
    #kk relations refractive index
    #integratedfn3(x) = (alpha_lb_2(x,E0_T_InAs_value,eps_inf,P,mhh,mlh,T,F,n_or_p) - alpha)/((x^2 - omega^2))
    #m_IB_p_value = eps_inf^(1.0/2.0) - 0.17 + c/pi*quadgk(integratedfn3,1.0,4.0/hbar_eV,rtol=1e-6)[1] # to make eps_inf at below bandgap frequencies (InAs ~0.06% error, InAs0.4SbP ~0.05% error)
    #constant refractive index 
    m_IB_p_value = eps_inf^(1/2)
    #dieletric form
    return (m_IB_p_value^2-m_IB_pp^2) + im*2*m_IB_p_value*m_IB_pp
end
precompile(epsIB,(Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64))



function epsIBEV(E_photon,N0,T,E0_T_InAs_value,eps_inf,P,mhh,mlh,F,n_or_p)
    #changing input eV to omega for photon number calculations
    omega = /(E_photon,hbar_eV)
    return epsIB(omega,N0,T,E0_T_InAs_value,eps_inf,P,mhh,mlh,F,n_or_p)
end
precompile(epsIBEV,(Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64))

end


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
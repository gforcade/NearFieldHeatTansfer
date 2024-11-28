__precompile__()
module ResponseModels
using FilmDataStructures,Optim
export eVtoOmega, eps_InAs_struct, eps_InAsntype, eps_InAsptype, hStep, cstRsp, sicRsp, siDsc, dptDsc, siRsp, prmMSi
# Collection of material permittivity models for use in heatLayers code.
# Conversion factors.
# Microns to electron volts.
const muEv = 1.239
# Radial frequency to electron volts.
const hBEv = 6.582e-16
# Temperature to electron volts.
const blzK = 8.6173e-5
# Electron volts to Joules.
const evJ = 1.602176565e-19
const eps_0 = 8.8542*^(10,-12)
const hbar = 1.05457*^(10,-34) #m2kg/s
const kb = 8.617333*^(10,-5) #eV K-1
const m_0 = 0.9109*^(10,-30)
const e = 1.6e-19 #C

@inline function cstRsp(cst::ComplexF64, enr::Float64)::ComplexF64

	return cst
end
precompile(cstRsp, (ComplexF64, Float64))
## Heaviside function for testing
@inline function hStep(arg::Float64)::Float64

	if arg >= 0.0

		return 1.0
	else
		return 0.0
	end
end
precompile(hStep, (Float64, ))
### Silicon carbide.
"""

	sicRsp(enr::Float64)::ComplexF64

Simplified permittivity response model for silicon carbide. 
Approximately valid for energies between 0.1 and 10 electron volts. 
"""
# Reproduction of the model reported in J. Appl. Phys. 106, 044306 (2009).
@inline function sicRsp(enr::Float64)::ComplexF64

	# Energy values (in electron volts) for the longitudinal and transverse optical phonons, 
	# and decay rate. 
	enrL = 0.12014342321784009
	enrT = 0.09831711568536007
	enrD = 0.0005901528146448003

	return 6.7 * (1.0 + /(enrL^2.0 - enrT^2.0, enrT^2.0 - enr^2.0 - im * enrD * enr))
end

precompile(sicRsp, (Float64, ))
### Silicon 
# Sample usage
# Boron p-type doping
# brnDpt = dptDsc(0.044, 1.0e18)
# Empty doping
# nllDpt = dptDsc(0.0, 0.0)
# Calculate silicon model parameters
# siModE = prmMSi(tmpLst[1], nllDpt, nllDpt)
# siRspE(enr) = siRsp(enr, siModE)
"""
Data structure holding computed model parameters for the permittivity response of silicon for 
given dopant properties at a fixed temperature. Call prmMSi(tmp::Float64, dnr::dptDsc, 
acp::dptDsc)::siDsc for construction. 

# Arguments
.tmp: temperature at which parameters where calculated. 
.dnrFac: strength of donor impurity Lorentzian.
.acpFac: strength of acceptor impurity Lorentzian.
.dnrDcy: decay rate (gamma) for donor impurity Lorentzian.
.acpDcy: decay rate (gamma) for acceptor impurity Lorentzian.
"""
struct siDsc

	tmp::Float64
	dnrFac::Float64
	acpFac::Float64
	dnrDcy::Float64
	acpDcy::Float64
end
# Include model code.
include("siMod.jl")
"""

	siRsp(enr::Float64, tmp::Float64, siVal::siDsc)::ComplexF64

Permittivity response model for doped silicon. See source file for additional details. 
"""
@inline function siRsp(enr::Float64, siVal::siDsc)::ComplexF64

	return siRspBck(enr, siVal.tmp) - /(siVal.dnrFac, enr^2 + im * enr * siVal.dnrDcy) - /(siVal.acpFac, enr^2 + im * enr * siVal.acpDcy)
end
precompile(siRsp, (Float64, siDsc))


using eps_InAs_ntype_v2
using eps_InAs_ptype
struct InAsDsc
    N0::Float64
    T::Float64
	mstar::Float64
	mstar_ptype::Float64
	mstar_ptype_hh::Float64
	mstar_ptype_lh::Float64
	gamma_ntype::Float64
	gamma_ptype::Float64
	E0_T_InAs_value::Float64
	E0_T_InAs_value_ptype::Float64
	eps_inf::Float64
	P::Float64
	F::Float64
end

@inline function eVtoOmega(energy) #energy in eV
    return(energy*evJ/hbar)
end
precompile(eVtoOmega, (Float64,))

@inline function eps_InAs_struct(N0,T)
	E0_T_InAs_value = E0_T_InAs(N0,T)
	E0_T_InAs_value_ptype = E0_T_InAs(N0,T)
	mstar= 0.024*m_0  
	del = 0.39  # eV, split-off band energy
	#mstar = m_star(N0,E0_T_InAs_value,m0_star,del)
	mstar_ptype_lh = 0.026*m_0
	mstar_ptype_hh = 0.36*m_0
	mstar_ptype = m_star_ptype(mstar_ptype_hh,mstar_ptype_lh)
	gamma_ntype = Gamma_ntype(N0,T,mstar)
	gamma_ptype = Gamma_ptype(N0,T,mstar_ptype)
	eps_inf = 11.6 #eps_inf parameter, same as in the gamma function
	P = 9.05e-8 #8.58e-8		#eV cm #comes from same place as InAsSbP params 
	#calculating the fermi energy
	Nc = 2.0*(3.0*E0_T_InAs_value*kb*T/(8.0*pi*P^2.0))^(3.0/2.0) #9.3301e16 #cm-3 #value from Sentaurus
	res = optimize(t -> fermi(t,1.0/(kb*T*np_factor(del,E0_T_InAs_value,mstar)),pi^(1/2)/2*N0/(Nc*10.0^6)),[E0_T_InAs_value],Newton()) #Nc*1e6 because: from cm-3 to m-3
	F = Optim.minimizer(res)[1]*kb*T + E0_T_InAs_value #fermi energy relative to valence band
    return InAsDsc(N0,T,mstar,mstar_ptype,mstar_ptype_hh,mstar_ptype_lh,gamma_ntype,gamma_ptype,E0_T_InAs_value,E0_T_InAs_value_ptype,eps_inf,P,F)
end
precompile(eps_InAs_struct,(Float64,Float64))

@inline function eps_InAsntype(E_photon::Float64,InAsDsc_struct::InAsDsc)::ComplexF64 #N_base in m^(-3)
     omega = eVtoOmega(E_photon)
     #o_angular = 5.38*10^14   #hbar*5.38*10^14/eV = 0.35eV
    #if(E_photon>InAsDsc_struct.E0_T_InAs_value)  #made just InAsDsc_struct be passed vs all the variables
	return epsIB(omega,InAsDsc_struct.N0,InAsDsc_struct.T,InAsDsc_struct.E0_T_InAs_value,InAsDsc_struct.eps_inf,InAsDsc_struct.P,InAsDsc_struct.mstar_ptype_hh,InAsDsc_struct.mstar_ptype_lh,InAsDsc_struct.F,0.0) + epsFCL(omega,InAsDsc_struct.mstar,InAsDsc_struct.gamma_ntype,InAsDsc_struct.N0,InAsDsc_struct.F-InAsDsc_struct.E0_T_InAs_value,InAsDsc_struct.E0_T_InAs_value,InAsDsc_struct.T,0) - InAsDsc_struct.eps_inf #eps_inf to avoid double counting
	#else
    #   return epsFCL(omega,InAsDsc_struct.mstar,InAsDsc_struct.gamma_ntype,InAsDsc_struct.N0)
    #end
end
precompile(eps_InAs_ntype_v2,(Float64,InAsDsc))


@inline function eps_InAsptype(E_photon::Float64,InAsDsc_struct::InAsDsc)::ComplexF64 #N_base in m^(-3)
	omega = eVtoOmega(E_photon)
	#o_angular = 5.38*10^14   #hbar*5.38*10^14/eV = 0.35eV
	#if(E_photon>InAsDsc_struct.E0_T_InAs_value_ptype)
	return epsIB(omega,InAsDsc_struct.N0,InAsDsc_struct.T,InAsDsc_struct.E0_T_InAs_value_ptype,InAsDsc_struct.eps_inf,InAsDsc_struct.P,InAsDsc_struct.mstar_ptype_hh,InAsDsc_struct.mstar_ptype_lh,InAsDsc_struct.F,1.0) + epsFCL(omega,InAsDsc_struct.mstar_ptype,InAsDsc_struct.gamma_ptype,InAsDsc_struct.N0,InAsDsc_struct.F-InAsDsc_struct.E0_T_InAs_value,InAsDsc_struct.E0_T_InAs_value,InAsDsc_struct.T,1) - InAsDsc_struct.eps_inf #eps_inf to avoid double counting
	#else
	#   return epsFCL_ptype(omega,InAsDsc_struct.mstar_ptype,InAsDsc_struct.gamma_ptype,InAsDsc_struct.N0)
    #end
end
precompile(eps_InAs_ptype,(Float64,InAsDsc))
end





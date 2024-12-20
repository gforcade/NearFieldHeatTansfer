# Model permittivity of doped silicon, reproducing the model used in Fu, C. J., 
# & Zhang, Z. M. (2006). Nanoscale radiation heat transfer for silicon at different doping levels. 
# (International Journal of Heat and Mass Transfer, 49(9-10)), and 1703-1718 American Journal of 
# Physics 44, 353 (1976); doi: 10.1119/1.10194.
## Dependencies
using Roots

# Conversion factors are loaded from call through filmUtilities.jl
"""

	blkFD(arg::Float64)::Float64

Blakemore approximation of the Fermi-Dirac integral, accurate to with two percent of true values. 
"""
function blkFD(arg::Float64)::Float64

	if arg < -5.0 

		return exp(arg)
	elseif arg < 0.7

		return exp(arg) - 0.3 * exp(2.0 * arg) + 0.06 * exp(3.0 * arg)
	elseif arg < 1.1

		return exp(arg) * /(1.0, 1.0 + 0.27 * exp(arg))
	elseif arg < 1.3

		return /(4.0 * ^(arg, 1.5), 3.0 * sqrt(pi)) * (1.0 + /(1.15, arg^2))
	else

		return /(4.0 * ^(arg, 1.5), 3.0 * sqrt(pi)) * (1.0 + /(pi^2, 8.0 * arg^2))
	end
end
precompile(blkFD, (Float64, ))


"""
	
	siRspBck(enr::Float64, tmp::Float64)::ComplexF64

Temperature dependent background permittivity of intrinsic silicon. Loss is assumed to be 
negligible in comparison to free carrier contributions. 
"""
# Reproduction of the model reported in Lee, B. J., Zhang, Z. M., Early, E. A., DeWitt, D. P., & 
# Tsai, B. K. (2005). Modeling radiative properties of silicon with coatings and comparison with 
# reflectance measurements. Journal of thermophysics and heat transfer, 19(4), 558-565. 
@inline function siRspBck(enr::Float64, tmp::Float64)::ComplexF64
	#Gavin: Changed T -> T - 273.15  #equation from reference says T is in celsius 
	#Timans
	kIB = SialphaIB(enr, tmp) * muEv * 1e-4 / (4.0 * pi * enr)   #multiplying by 1e-4 to convert from um to cm ## verified
	#J-M
	#return ^(sqrt(4.565 + /(97.3, ^(3.648,2) - ^(/(enr * 1.24, muEv), 2))) + (-1.864 * ^(10, -4) + 5.394 * /(^(10, -3), ^(3.648,2) - ^(/(enr * 1.24, muEv), 2))) * (tmp - 273.15) + im*kIB, 2) #- /(^(hBEv * 1.0e13, 2), enr * (enr + im * 1.0e10 * hBEv)) 
	#Li #verified
	lm = muEv / enr #lmabda in um
	eps_r = 11.4445 + 2.7739e-4 * tmp + 1.7050e-6 * tmp^2.0 - 8.1347e-10 * tmp^3.0
	g = 0.8948 + 4.3977e-4 * tmp + 7.3835e-8 * tmp^2.0 
	n = exp(-3.0*(-0.071 + 1.887e-6 * tmp + 1.934e-9 * tmp^2.0 - 4.544e-13 * tmp^3.0))
	return ^( sqrt(eps_r + g*n / lm^2.0) + im*kIB,2.0)
end
precompile(siRspBck, (Float64, Float64))

"""
	SikIB(enr::Float64,tmp::Float64)::Float64

	Temperature dependent background permittivity of intrinsic silicon.  From Timans, Emissivity of silicon at elevated temperatures (1993)
"""
@inline function SialphaIB(enr::Float64, tmp::Float64)::Float64


	blzk = 8.6173e-5
	return 1.0/enr * (F1(enr-siEnrGap(tmp)+blzk*212.0)/(exp(212.0/tmp) - 1.0) + F1(enr-siEnrGap(tmp)-blzk*212.0)/(1.0 - exp(-212.0/tmp)) + F2(enr-siEnrGap(tmp)+blzk*670.0)/(exp(670.0/tmp) - 1.0) + F2(enr-siEnrGap(tmp)-blzk*670.0)/(1.0 - exp(-670.0/tmp)) ) + Sialpha3a(enr,tmp) + Sialpha4a(enr,tmp)
end
precompile(SialphaIB, (Float64, Float64 ))

@inline function F1(x::Float64)::Float64

	if x > 0.0055
		
		return 0.504*sqrt(x) + 392.0*(x - 0.0055)^2.0
	elseif x > 0.0
		 
		return 0.504*sqrt(x)
	else
		return 0.0
	end
end
precompile(F1, (Float64, ))

@inline function F2(x::Float64)::Float64

	if x > 0.0055
		
		return 18.08*sqrt(x) + 5760.0*(x - 0.0055)^2.0
	elseif x > 0.0
		 
		return 18.08*sqrt(x)
	else
		return 0.0
	end
end
precompile(F2, (Float64, ))

@inline function Sialpha3a(enr::Float64, tmp::Float64)::Float64
	#alpha 3a

	return Sialphaia(enr,tmp,536.0,1050.0)
end
precompile(Sialpha3a, (Float64, Float64 ))

@inline function Sialpha4a(enr::Float64, tmp::Float64)::Float64
	#alpha 4a

	return Sialphaia(enr,tmp,988.0,1420.0)
end
precompile(Sialpha4a, (Float64, Float64 ))

@inline function Sialphaia(enr::Float64, tmp::Float64, coef::Float64, Tcoeff::Float64)::Float64
	#alpha for 3a and 4a


	if enr < siEnrGap(tmp) - 8.6173e-5 * Tcoeff
		return 0.0
	else
		return coef*(enr - siEnrGap(tmp) + 8.6173e-5*Tcoeff)^2.0 / (enr*(exp(Tcoeff/tmp) - 1.0))
	end
end
precompile(Sialphaia, (Float64, Float64, Float64, Float64 ))
"""

	siEnrGap(tmp::Float64)::Float64 

Temperature dependent energy gap for intrinsic silicon in electron volts.
V. Alex, S. Finkbeiner, and J. Weber, “Temperature dependence of the indirect energy gap in crystalline silicon ARTICLES YOU MAY BE INTERESTED IN,” J. Appl. Phys., vol. 79, p. 6943, 1996, doi: 10.1063/1.362447
"""
@inline function siEnrGap(tmp::Float64)::Float64 

	return 1.1692 - /(0.00049 * tmp^2.0, tmp + 655.0)
end
precompile(siEnrGap, (Float64, ))
"""

	itrnV(tmp::Float64)::Float64

Number of intrinsic valence band carriers in units of cm^{-3} as a function of temperature in 
Kelvin. 
"""
@inline function itrnV(tmp::Float64)::Float64

	return 2.66 * ^(10.0, 19.0) * ^(/(tmp, 300), 1.5)
end
precompile(itrnV, (Float64, ))
"""

	itrnC(tmp::Float64)::Float64

Number of intrinsic conduction band carriers in units of cm^{-3} as a function of temperature in 
Kelvin. 
"""
@inline function itrnC(tmp::Float64)::Float64 

	return 2.86 * ^(10.0, 19.0) * ^(/(tmp, 300), 1.5)
end
precompile(itrnC, (Float64, ))
"""

	SiMobn(dnr::dptDsc)::Float64 

Electron mobility vs doping concentration at room temperature. Doping in cm^-3. Output in cm^2/(V s)
From: Basu et al. "Infrared radiative properties of heavily doped silicon at room temperature", 2010
"""
@inline function SiMobn(eCon::Float64)::Float64 

	return 68.5 + (1414.0 - 68.5) / (1.0 + (eCon/9.2e16)^0.711) - 56.5 / (1.0 + (3.41e20 / eCon)^1.98) 
end
precompile(SiMobn, (Float64, ))
"""

	SiMobp(acp::dptDsc)::Float64 


Hole mobility vs doping concentration at room temperature. Doping in cm^-3. Output in cm^2/(V s)
From: Basu et al. "Infrared radiative properties of heavily doped silicon at room temperature", 2010
"""
@inline function SiMobp(hCon::Float64)::Float64 

	return 44.9 * exp(-9.23e16/hCon) + 470.0 / (1.0 + (hCon/9.2e16)^0.719) - 29.0 / (1.0 + (3.41e20 / hCon)^2.0) 
end
precompile(SiMobp, (Float64, ))
"""

	function dnrCcn(enrFrm::Float64, tmp::Float64, acp::dptDsc)::Float64	

Acceptor concentration function. 
"""
@inline function acpCcn(enrFrm::Float64, tmp::Float64, acp::dptDsc)::Float64

	return acp.ccn * (1.0 - /(4.0, 4.0  + exp(/(enrFrm - acp.enr, blzK * tmp))))
end
precompile(acpCcn, (Float64, Float64, dptDsc))
"""

	function dnrCcn(enrFrm::Float64, tmp::Float64, dnr::dptDsc)::Float64

Donor concentration function. 
"""
@inline function dnrCcn(enrFrm::Float64, tmp::Float64, dnr::dptDsc)::Float64

	return dnr.ccn * (1.0 - /(2.0, 2.0  + exp(/(siEnrGap(tmp) - dnr.enr - enrFrm, blzK * tmp))))
end
precompile(dnrCcn, (Float64, Float64, dptDsc))
"""
	
	frmZero(enrFrm::Float64, tmp::Float64, dnr::dptVal, acp::dptVal)::Float64

Net charge concentration as a function of supposed Fermi energy and temperature. 
Used as an implicit equation for determining the correct Fermi energy. 
"""
@inline function frmZero(enrFrm::Float64, tmp::Float64, dnr::dptDsc, acp::dptDsc)::Float64

	#return itrnV(tmp) * blkFD(-/(enrFrm, blzK * tmp)) + /(dnr.ccn, 1.0 + 2.0 * exp(/(enrFrm - dnr.enr, blzK * tmp))) - itrnC(tmp) * blkFD(/(enrFrm - siEnrGap(tmp), blzK * tmp)) - /(acp.ccn, 1.0 + 4.0 * exp(/(acp.enr - enrFrm, blzK * tmp)))
	return itrnV(tmp) * blkFD(-/(enrFrm, blzK * tmp)) + dnrCcn(enrFrm,tmp,dnr) - itrnC(tmp) * blkFD(/(enrFrm - siEnrGap(tmp), blzK * tmp)) - acpCcn(enrFrm,tmp,acp)

end
precompile(frmZero, (Float64, Float64, dptDsc, dptDsc))
"""

	function Ci(enrFrm::Float64, tmp::Float64, dnr::dptDsc)::Float64
Basu, Lee, Zhang. "Infrared radiative properties of heavily doped silicon at room temperature", 2010.
Ratio of concentration of ionized dopant atoms to the total doping concentration
Used this model in Mittapaly et al. 2021. Provided excellent fitting.
"""
@inline function Ci( tmp::Float64, N::Float64, Dopant_Type::Bool)::Float64
	#N is dopant concentration
	#tmp is temperature of material
	if Dopant_Type == true
		#n-type doping
		A = 0.0824 * (tmp / 300.0)^(-1.622)
		N0 = 1.6e18 * (tmp / 300.0)^(0.7267)
		if N < N0
			B = 0.4722 * (tmp / 300.0)^(0.0652)
		else
			B = 1.23 - 0.3162 * (tmp / 300.0)
		end
	else
		#p-type doping
		A = 0.2364 * (tmp / 300.0)^(-1.474)
		N0 = 1.577e18 * (tmp / 300.0)^(0.46)
		if N < N0
			B = 0.433 * (tmp / 300.0)^(0.2213)
		else
			B = 1.268 - 0.338 * (tmp / 300.0)
		end
	end
	return 1.0 - A*exp(-(B*log(N/N0))^2.0)
end
precompile(dnrCcn, (Float64, Float64, Bool))
"""

	prmMSi(tmp::Float64, dnr::dptDsc, acp::dptDsc)::siDsc

Computation of model parameters for doped silicon, see preamble of source file for details. 
"""
function prmMSi(tmp::Float64, dnr::dptDsc, acp::dptDsc)::siDsc

	# Self-consistent Fermi level.
	frmFnc(flvl) = frmZero(flvl, tmp, dnr, acp)					#verified
	enrFrm = find_zero(frmFnc, /(siEnrGap(tmp), 2.0))			#verified
	
	
	# assumes model from Gaylord
	# Number of ionized acceptor and donor dopants (carriers), units cm^-3. Does not include free carriers of undoped material
	#cncDnr = dnrCcn(enrFrm, tmp, dnr)							#verified
	#cncAcp = acpCcn(enrFrm, tmp, acp)							#verified
	"""
	# Number of free electrons and holes (including from ionized dopants) Better than before as it now can accurately count carriers for low doping
	cncDnr = itrnC(tmp) * blkFD(/(enrFrm - siEnrGap(tmp), blzK * tmp))
	cncAcp = itrnV(tmp) * blkFD(-/(enrFrm, blzK * tmp))
	"""

	## Assumes Basu 2010 model
	#Number of free electrons and holes
	Nth2 = itrnC(tmp) * itrnV(tmp) * exp(-/(siEnrGap(tmp), blzK * tmp))
	if acp.ccn == 0.0
		cncDnr = Ci(tmp,dnr.ccn,true) * dnr.ccn 
		cncDnr = 0.5 * (cncDnr + sqrt(cncDnr^2.0 + 4.0*Nth2))
		cncAcp = Nth2 / cncDnr
	else
		cncAcp = Ci(tmp,acp.ccn,false) * acp.ccn 
		cncAcp = 0.5 * (cncAcp + sqrt(cncAcp^2.0 + 4.0*Nth2))
		cncDnr = Nth2 / cncAcp
	end
	
	# Electron and hole effective masses. #verified
	mEffE = 0.27 
	mEffH = 0.37

	"""
	## Fu and Zhang model
	# Room temperature total electron scattering time.
	stETR = (141.0 + /(19.5, 1.0 + ^(/(dnr.ccn, 1.3 * ^(10.0, 17.0)), 0.91))) * 1.0e-15		#verified
	# Room temperature total hole scattering time.
	stHTR = (10.0 + /(94.0, 1.0 + ^(/(acp.ccn, 1.9 * ^(10.0, 17.0)), 0.76))) * 1.0e-15		#verified
	"""

	## Basu 2010 model
	# Room temperature total electron scattering time.
	stETR = ^(e / (mEffE*m_0 * SiMobn(dnr.ccn)/10000.0), -1.0)		#verified
	# Room temperature total hole scattering time.
	stHTR = ^(e / (mEffH*m_0 * SiMobp(acp.ccn)/10000.0),-1.0)		#verified
	
	# Room temperature electron-lattice scattering time.
	stELR = 2.23 * 1.0e-13											#verified
	# Room temperature hole-lattice scattering time.
	stHLR = 1.06 * 1.0e-13											#verified
	# Room temperature electron-impurity scattering time. 
	stEIR = /(stETR * stELR, stELR - stETR) 						#verified
	# Room temperature hole-impurity scattering time. 
	stHIR = /(stHTR * stHLR, stHLR - stHTR) 						#verified
    
	# Include temperature dependence scattering times.
	# Electrons 
	stEI = stEIR * ^(/(tmp, 300.0), 1.5)								#verified
	stEL = stELR * ^(/(tmp, 300.0), -3.8)								#verified
	# Holes
	stHI = stHIR * ^(/(tmp, 300.0), 1.5)								#verified
	stHL = stHLR * ^(/(tmp, 300.0), -3.6)								#verified

	# Total donor and acceptor decay rates.  in eV
	dnrDcy = hBEv * (/(1.0, stEI) + /(1.0, stEL)) 
	acpDcy = hBEv * (/(1.0, stHI) + /(1.0, stHL)) 


	# Conversion factors for donor and acceptor contributions. units eV^2
	dnrFac = /(cncDnr * hBEv^2 * ^(10.0, 9.0), mEffE * 0.5199895 * 55.26349406 * 0.01112265)  #verified
	acpFac = /(cncAcp * hBEv^2 * ^(10.0, 9.0), mEffH * 0.5199895 * 55.26349406 * 0.01112265) #verified

	return siDsc(tmp, dnrFac, acpFac, dnrDcy, acpDcy)
end
precompile(prmMSi, (Float64, dptDsc, dptDsc))
__precompile__()
# See readme for module description. 
module FilmUtilities
using Base.Threads, Cubature, ProgressMeter
## Load data structures. 
using FilmDataStructures
export plcFnc, tfrIntEv!, heatTfr!
## Constants
# Relative tolerance for cubature. 1.0e-4, 1.0e-3
const cubRelTol = 1.0e-4
const cubAbsTolCoeff = 1.0e-3
# Cutoff for unbounded integrals tan and atan variable transforms.
# tan(pi/2 - intDltInf) scaling. 
# Use number smaller than 0.5e-4 for sub nanometer gaps. 1.0e-3. 0.5e-4, no change. 1e-2, low energy lower.
const intDltInf = 1e-5
# Initial guess for width of ``intermediate wave'' behavior. 1.0e-3
const lcWdt = 1.0e-3
# Desired relative magnitude for defining width of ``intermediate wave'' behavior. 5.0e-2
const lcRel = 5.0e-2
# Ratio of /(real(rsp(ener)), imag(rsp(ener))) before intermediate wave point is created. 1.0e3
const iTrnV = 1.0e3
# Near light line normalization, effectively number of wavelengths before decay.50.0. 1e3, slightly better
const mxSz = 50.0
# Near light regularization, imposed as imaginary part of perpendicular wave vector. 1.0e-2
const lcReg = 1.0e-2
# Conversion factors.
# Microns to electron volts.
const muEv = 1.239
# Radial frequency to electron volts.
const hBEv = 6.582e-16
# Temperature to electron volts.
const blzK = 8.6173e-5
# Electron volts to Joules.
const evJ = 1.602176565e-19
# Repository functions no longer being used, but possibly helpful for testing. 
# include("filmUtilitiesRepo.jl")
"""

	geoSrs(intTrm::ComplexF64, geoFct::ComplexF64)::ComplexF64

Analytic extension of the geometric series for a given initial term and geometric factor.
"""
@inline function geoSrs(intTrm::ComplexF64, geoFct::ComplexF64)::ComplexF64

	return /(intTrm, 1.0 - geoFct)
end
precompile(geoSrs, (ComplexF64, ComplexF64))
"""

	zSqrt(val::ComplexF64)::ComplexF64

Square root defined with cut line along real axis. 
"""
@inline function zSqrt(val::ComplexF64)::ComplexF64

	if angle(val) >= 0.0

		return sqrt(val)
	else

		return sqrt(val) * exp(im * pi)
	end
end
precompile(zSqrt, (ComplexF64, ))
"""

	wvcPrp(rsp::ComplexF64, lTck::Float64, wvc::Float64)::ComplexF64

Normalized perpendicular wave, incorporating an effective normalization.
"""
@inline function wvcPrp(rsp::ComplexF64, lTck::Float64, wvc::Float64)::ComplexF64

	return zSqrt(rsp - wvc^2) + im * lcReg * ^(/(sqrt(abs(zSqrt(rsp - wvc^2))^2 + wvc^2) * lTck, mxSz), 4)
end
precompile(wvcPrp, (ComplexF64, Float64, Float64))
"""

	plcFnc(enr::Float64, tmp::Float64)::Float64

Planck function for energy in electron volts and temperature in Kelvin.
"""
@inline function plcFnc(enr::Float64, tmp::Float64)::Float64

	return /(enr, exp(/(enr, blzK * tmp)) - 1.0)
end
precompile(plcFnc, (Float64, Float64))
"""
	flxPfc(lVar::lyrDsc, lPair::Union{Array{Int64,1},SubArray{Int64,1}}, enr::Float64)::Float64

Energy prefactor, W eV^-1 cm^-2, for fluctuating electric field intensity integrand. 
"""
@inline function flxPfc(lVar::lyrDsc, lPair::Union{Array{Int64,1},SubArray{Int64,1}}, enr::Float64)::Float64

	#return /(evJ * enr^2 * pi^3, 2.0^2 * muEv^2 * hBEv * ^(10.0, -8)) * (plcFnc(enr, lVar.tmpLst[lPair[1]]) - plcFnc(enr, lVar.tmpLst[lPair[2]])) ###error in units
	return /(evJ * enr^2 * 2.0, muEv^2 * hBEv * ^(10.0, -8)) * (plcFnc(enr, lVar.tmpLst[lPair[1]]) - plcFnc(enr, lVar.tmpLst[lPair[2]]))
end
precompile(flxPfc, (lyrDsc, Union{Array{Int64,1},SubArray{Int64,1}}, Float64))	
"""

	prpLyr(dst::Float64, rsp::ComplexF64, wvc::Float64)::ComplexF64

Calculation of phase accumulation under propagation, distance assumed to be in units of wavelength.
"""
@inline function prpLyr(dst::Float64, rsp::ComplexF64, wvc::Float64)::ComplexF64

	return exp(2.0 * pi * im * wvcPrp(rsp, dst, wvc) * dst)
end
precompile(prpLyr, (Float64, ComplexF64, Float64))
"""
	lyrTck(lVar::lyrDsc, lNum::Int)::Float64

Layer thickness in microns. 
"""
@inline function lyrTck(lVar::lyrDsc, lNum::Int)::Float64

	if lNum == 1 || lNum == length(lVar.tmpLst)

		return 0.0
	else

		return lVar.bdrLoc[lNum] - lVar.bdrLoc[lNum - 1]
	end
end
precompile(lyrTck, (lyrDsc, Int))
"""

	lyrTckRel(lVar::lyrDsc, lNum::Int, enr::Float64)::Float64

Compute layer thickness in units of relative wavelength.
"""
@inline function lyrTckRel(lVar::lyrDsc, lNum::Int, enr::Float64)::Float64

	return /(enr * lyrTck(lVar, lNum), muEv)
end
precompile(lyrTckRel, (lyrDsc, Int, Float64))
"""

	strTck(lVar::lyrDsc, enr::Float64)::Float64

Compute total structure thickness in units of relative wavelength.
"""
@inline function strTck(lVar::lyrDsc, enr::Float64)::Float64

	return /(enr * (lVar.bdrLoc[end] - lVar.bdrLoc[1]), muEv)
end
precompile(strTck, (lyrDsc, Float64))
"""

	xfrBdr!(kPrp::NTuple{2,ComplexF64}, rsp::NTuple{2,ComplexF64}, bdrArr::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing

Isotopic electric reflection and transmission coefficients for two half spaces. 
"""
# First entries in the permittivity (rps) and wave vector (kPrp) tuples are consider to be 
# incident. 
# The first two entries returned are for 'p' polarized waves, the next two for 's' 
# polarized waves. 
# Odd numbers are transmission, even modes are reflection.
@inline function xfrBdr!(kPrp::NTuple{2,ComplexF64}, rsp::NTuple{2,ComplexF64}, bdrArr::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing
	# Fresnel coefficients.
	copyto!(bdrArr, [/(2.0 * zSqrt(rsp[1] * rsp[2]) * kPrp[1], rsp[1] * kPrp[2] + rsp[2] * kPrp[1]), /(rsp[1] * kPrp[2] - rsp[2] * kPrp[1], rsp[1] * kPrp[2] + rsp[2] * kPrp[1]), /(2.0 * kPrp[1], kPrp[2] + kPrp[1]), /(kPrp[2] - kPrp[1], kPrp[2] + kPrp[1])])

	 return nothing
end
precompile(xfrBdr!, (NTuple{2,ComplexF64}, NTuple{2,ComplexF64}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}))
"""

	bdrFill!(lVar::lyrDsc, enr::Float64, wvc::Float64, mode::Int, bdrArr::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}})::Nothing

Fills bdrArr with individual boundary transfer coefficients, i.e. as though every boundary were 
between half-spaces. 
"""
# If mode == 1 smaller indices are treated as initial. 
# If mode == 2 larger indices are treated as initial.
# Convention for a given column of bdrArr as described for xfrBdr. 
@inline function bdrFill!(lVar::lyrDsc, enr::Float64, wvc::Float64, mode::Int, bdrArr::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}})::Nothing
	# Smaller indices are treated as initial. 
	if mode == 1

		for bdr = 1 : length(lVar.bdrLoc)

			@views xfrBdr!((wvcPrp(lVar.rspPrf[bdr](enr), lyrTckRel(lVar, bdr, enr), wvc), wvcPrp(lVar.rspPrf[bdr + 1](enr), lyrTckRel(lVar, bdr + 1, enr), wvc)), (lVar.rspPrf[bdr](enr), lVar.rspPrf[bdr + 1](enr)), bdrArr[:,bdr])
		end
	# Larger indices are treated as initial. 
	elseif mode == 2

		for bdr = 1 : length(lVar.bdrLoc)
			
			@views xfrBdr!((wvcPrp(lVar.rspPrf[bdr + 1](enr), lyrTckRel(lVar, bdr + 1, enr), wvc), wvcPrp(lVar.rspPrf[bdr](enr), lyrTckRel(lVar, bdr, enr), wvc)), (lVar.rspPrf[bdr + 1](enr), lVar.rspPrf[bdr](enr)), bdrArr[:,bdr])
		end
	else

		error("Unrecognized fill mode.")	
	end

	return nothing 
end
precompile(bdrFill!, (lyrDsc, Float64, Float64, Int, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}))
"""
	
	trnFlmExp!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, trnArrHL::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, trnArrHR::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, refArrHL::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, refArrHR::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, trnArr::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing

Updates transmission coefficients under expansion of nested reflection coefficient. 
See associated notes for additional details. 
"""
# If mode == 1, expansion is ``to the right'', i.e. expanding to growing layers numbers. 
# If mode == 2, expansion is ``to the left'', i.e. expanding to shrinking layers numbers.
# First entry of trnArr should be numerator coefficient. 
# The next two entries are the reflection and additional coefficients for the denominator. 
# See notes for additional details. 
@inline function trnFlmExp!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, trnArrHL::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, trnArrHR::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, refArrHL::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, refArrHR::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, trnArr::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing

	if mode == 1
	
		copyto!(trnArr, [trnArrHL[lNum - 1] * prpLyr(lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * trnArr[1], trnArr[2] + refArrHL[lNum - 1] * trnArr[3], prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * (trnArr[3] - refArrHR[lNum - 1] * trnArr[2])])
	elseif mode == 2
		
		copyto!(trnArr, [trnArrHR[lNum] * prpLyr(lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * trnArr[1], trnArr[2] + refArrHR[lNum] * trnArr[3], prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum ](enr), wvc) * (trnArr[3] - refArrHL[lNum] * trnArr[2])])
	else 

		error("Unrecognized increment mode.")	
	end

	return nothing
end
precompile(trnFlmExp!, (lyrDsc, Float64, Float64, Int, Int, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}))
"""
	
	xfrFlmNst!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, refArrHL::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, refArrHR::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, xfrArr::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing

Updates transfer coefficients by nesting current values, i.e. inclusion of an additional layer. 
See associated notes for additional details. 
"""
# If mode == 1, nesting is ``to the right'', i.e. adding growing layer numbers. 
# If mode == 2, nesting is ``to the left'', i.e. adding shrinking layer numbers.
# xfrArr should hold the expansion coefficients, to the current level, of the numerator as the 
# its first two entries, additional term followed by reflection term. 
# The next two entries follow the same convention, but for the denominator. 
# See notes for additional details. 
@inline function xfrFlmNst!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, refArrHL::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, refArrHR::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, xfrArr::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing

	if mode == 1
	
		copyto!(xfrArr, [refArrHR[lNum] * xfrArr[3] + prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * xfrArr[1], refArrHR[lNum] * xfrArr[4] + prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * xfrArr[2], xfrArr[3] - prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * refArrHL[lNum] * xfrArr[1], xfrArr[4] - prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * refArrHL[lNum] * xfrArr[2]])
	elseif mode == 2
		
		copyto!(xfrArr, [refArrHL[lNum - 1] * xfrArr[3] + prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * xfrArr[1], refArrHL[lNum - 1] * xfrArr[4] + prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * xfrArr[2], xfrArr[3] - prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * refArrHR[lNum - 1] * xfrArr[1], xfrArr[4] - prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * refArrHR[lNum - 1] * xfrArr[2]])
	else 

		error("Unrecognized increment mode.")	
	end

	return nothing
end
precompile(xfrFlmNst!, (lyrDsc, Float64, Float64, Int, Int, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}))
"""

	relFlm!(lVar::lyrDsc, enr::Float64, wvc::Float64, lPair::Union{Array{Int64,1},SubArray{Int64,1}}, imdCff::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, xfrCff::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing

Relative field amplitude coefficients from transfer coefficients. 
"""
# imdCff[:, lyr] = [trnCf; refLeftCellLeft; refLeftCellRight; refRightCellLeft; 
# refRightCellRight], with p-pol followed by s-pol in each. 
function relFlm!(lVar::lyrDsc, enr::Float64, wvc::Float64, lPair::Union{Array{Int64,1},SubArray{Int64,1}}, imdCff::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, xfrCff::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing
	# Layer phase propagation. 
	if lPair[1] > 1

		lLyrPhz = prpLyr(lyrTckRel(lVar, lPair[1], enr), lVar.rspPrf[lPair[1]](enr), wvc)
	else

		lLyrPhz = 0.0
	end

	if lPair[2] < length(lVar.tmpLst)

		rLyrPhz = prpLyr(lyrTckRel(lVar, lPair[2], enr), lVar.rspPrf[lPair[2]](enr), wvc)
	else

		rLyrPhz = 0.0
	end
	# Field dressing, absolute value used for stability, consistent with later expressions. 
	pDrs = abs(geoSrs(1.0 + im * 0.0, imdCff[3] * lLyrPhz * imdCff[5] * lLyrPhz)) * abs(geoSrs(1.0 + im * 0.0, imdCff[7] * rLyrPhz * imdCff[9] * rLyrPhz))
	sDrs = abs(geoSrs(1.0 + im * 0.0, imdCff[4] * lLyrPhz * imdCff[6] * lLyrPhz)) * abs(geoSrs(1.0 + im * 0.0, imdCff[8] * rLyrPhz * imdCff[10] * rLyrPhz))
	# Propagation term that should be included in left moving waves is shifted to the 
	# flux integrand to avoid overflow. 
	copyto!(xfrCff, [pDrs * imdCff[1], imdCff[9] * imdCff[1], sDrs * imdCff[2], sDrs * imdCff[10] * imdCff[2]])

	return nothing
end
precompile(relFlm!, (lyrDsc, Float64, Float64, Union{Array{Int64,1},SubArray{Int64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}))
"""
	
	tfrFlm!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, wvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}})::Nothing

Calculates field transfer coefficients for a list of ordered pairs of layers. 
"""
# bdrArrL will contain interface transfer coefficients with smaller layer numbers treated as 
# initial. 
# bdrArrR will contain interface transfer coefficients with larger layer numbers treated as 
# initial. 
function tfrFlm!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, wvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}})::Nothing
	# Calculate infinite interface coefficients.
	# Left incident.
	bdrFill!(lVar, enr, wvc, 1, bdrArrL)
	# Right incident.
	bdrFill!(lVar, enr, wvc, 2, bdrArrR)
	# Double check that pairs in list are correctly ordered. 
	for pr = 1:size(lPairs)[2]

		if lPairs[2, pr] <= lPairs[1, pr]
			
			error("Current program convention requires ordered layers.")
			return nothing
		end
	end
	# Unique, sorted, source and target layer lists.
	srcLyrs = sort!(unique(lPairs[1,:]))
	trgLyrs = sort!(unique(lPairs[2,:]))
	# Running indices for unique source and target layers.
	srcUInd = 1
	trgUInd = 1
	# Number of layers.
	numLyrs = length(lVar.tmpLst)
	# Reflection and transmission coefficient views. 
	@views refLP = bdrArrL[2, :]
	@views refLS = bdrArrL[4, :]
	@views refRP = bdrArrR[2, :]
	@views refRS = bdrArrR[4, :]
	@views trnLP = bdrArrL[1, :]
	@views trnLS = bdrArrL[3, :]
	@views trnRP = bdrArrR[1, :]
	@views trnRS = bdrArrR[3, :]
	# Temporary transmission coefficients.
	trnTmpP = [1.0 + im * 0.0, 1.0 + im * 0.0, 0.0 + im * 0.0]
	trnTmpS = [1.0 + im * 0.0, 1.0 + im * 0.0, 0.0 + im * 0.0]
	### Calculate layer reflection and transmission coefficients.
	## Right incident reflection coefficients. 
	# Safeguard against fictitious reflections from infinite boundary.
	if srcLyrs[1] == 1

		for lyrPInd in findall(x -> x == 1, lPairs[1, :])

			imdCff[3:4, lyrPInd] .= [0.0 + im * 0.0, 0.0 + im * 0.0]
		end
		# Move to next set of source layers. 
		if length(srcLyrs) > srcUInd
			
			srcUInd += 1
		end
	end 
	# Seed nesting construction.
	# P-polarized.
	xfrCffTmpP = [0.0 + im * 0.0, 1.0 + im * 0.0, 1.0 + im * 0.0, 0.0 + im * 0.0]
	# S-polarized.
	xfrCffTmpS = [0.0 + im * 0.0, 1.0 + im * 0.0, 1.0 + im * 0.0, 0.0 + im * 0.0]
	# Set index to begin nesting.
	ind = 2
	
	while ind < trgLyrs[end]
		# Save right incident reflection coefficients for source cells.
		if srcLyrs[srcUInd] == ind
			
			for lyrPInd in findall(x -> x == ind, lPairs[1, :])
				# P-polarized. 
				imdCff[3, lyrPInd] = /(xfrCffTmpP[1] + xfrCffTmpP[2] * refRP[1], xfrCffTmpP[3] + xfrCffTmpP[4] * refRP[1])
				# S-polarized. 
				imdCff[4, lyrPInd] = /(xfrCffTmpS[1] + xfrCffTmpS[2] * refRS[1], xfrCffTmpS[3] + xfrCffTmpS[4] * refRS[1])
			end
			# Move to next set of source layers. 
			if length(srcLyrs) > srcUInd
				
				srcUInd += 1	
			end
		end
		# Save right incident reflection coefficients for target cells.
		if length(trgLyrs) >= trgUInd && trgLyrs[trgUInd] == ind
			
			for lyrPInd in findall(x -> x == ind, lPairs[2, :])
				# P-polarized
				imdCff[7, lyrPInd] = /(xfrCffTmpP[1] + xfrCffTmpP[2] * refRP[1], xfrCffTmpP[3] + xfrCffTmpP[4] * refRP[1])
				# S-polarized
				imdCff[8, lyrPInd] = /(xfrCffTmpS[1] + xfrCffTmpS[2] * refRS[1], xfrCffTmpS[3] + xfrCffTmpS[4] * refRS[1])
			end
			# Move to next set of target layers. 
			if length(trgLyrs) > trgUInd
				
				trgUInd += 1	
			end
		end
		# Increment (nest) reflection field coefficients.
		# P-polarized.
		xfrFlmNst!(lVar, enr, wvc, ind, 1, refLP, refRP, xfrCffTmpP)
		# S-polarized.
		xfrFlmNst!(lVar, enr, wvc, ind, 1, refLS, refRS, xfrCffTmpS)
		ind = ind + 1
	end
	# Save right incident reflection coefficients for largest indexed target cells.
	for lyrPInd in findall(x -> x == trgLyrs[end], lPairs[2, :])
		# P-polarized.
		imdCff[7, lyrPInd] = /(xfrCffTmpP[1] + xfrCffTmpP[2] * refRP[1], xfrCffTmpP[3] + xfrCffTmpP[4] * refRP[1])
		# S-polarized.
		imdCff[8, lyrPInd] = /(xfrCffTmpS[1] + xfrCffTmpS[2] * refRS[1], xfrCffTmpS[3] + xfrCffTmpS[4] * refRS[1])
	end
	## Left incident reflection coefficients.
	# Safeguard against fictitious reflections from infinite boundary.
	if trgLyrs[end] == numLyrs

		for lyrPInd in findall(x -> x == numLyrs, lPairs[2, :])

			imdCff[9:10, lyrPInd] .= [0.0 + im * 0.0, 0.0 + im * 0.0]
		end
		# Move to next set of target layers.
		if trgUInd > 1
			
			trgUInd -= 1
		end
	end
	# Reseed nesting construction. 
	# P-polarized.
	xfrCffTmpP .= [0.0 + im * 0.0, 1.0 + im * 0.0, 1.0 + im * 0.0, 0.0 + im * 0.0]
	# S-polarized.
	xfrCffTmpS .= [0.0 + im * 0.0, 1.0 + im * 0.0, 1.0 + im * 0.0, 0.0 + im * 0.0]
	# Set index to begin nesting.
	ind = numLyrs - 1

	while ind > srcLyrs[1]
		# Save left incident reflection coefficients for source cells.
		if srcLyrs[srcUInd] == ind
			
			for lyrPInd in findall(x -> x == ind, lPairs[1, :])
				# P-polarized
				imdCff[5, lyrPInd] = /(xfrCffTmpP[1] + xfrCffTmpP[2] * refLP[end], xfrCffTmpP[3] + xfrCffTmpP[4] * refLP[end])
				# S-polarized
				imdCff[6, lyrPInd] = /(xfrCffTmpS[1] + xfrCffTmpS[2] * refLS[end], xfrCffTmpS[3] + xfrCffTmpS[4] * refLS[end])
			end
			# Move to next set of source layers. 
			if srcUInd > 1

				srcUInd -= 1	
			end
		end
		# Save left incident reflection coefficients for target cells.
		if trgLyrs[trgUInd] == ind
			
			for lyrPInd in findall(x -> x == ind, lPairs[2, :])
				# P-polarized
				imdCff[9, lyrPInd] = /(xfrCffTmpP[1] + xfrCffTmpP[2] * refLP[end], xfrCffTmpP[3] + xfrCffTmpP[4] * refLP[end])
				# S-polarized
				imdCff[10, lyrPInd] = /(xfrCffTmpS[1] + xfrCffTmpS[2] * refLS[end], xfrCffTmpS[3] + xfrCffTmpS[4] * refLS[end])
			end
			# Move to next set of target layers. 
			if trgUInd > 1
				
				trgUInd -= 1	
			end
		end
		# Increment reflection field coefficients.
		# P-polarized
		xfrFlmNst!(lVar, enr, wvc, ind, 2, refLP, refRP, xfrCffTmpP)
		# S-polarized
		xfrFlmNst!(lVar, enr, wvc, ind, 2, refLS, refRS, xfrCffTmpS)
		ind = ind - 1
	end
	# Save left incident reflection coefficients for smallest indexed source cells.
	for lyrPInd in findall(x -> x == srcLyrs[1], lPairs[1, :])
		# P-polarized
		imdCff[5, lyrPInd] = /(xfrCffTmpP[1] + xfrCffTmpP[2] * refLP[end], xfrCffTmpP[3] + xfrCffTmpP[4] * refLP[end])
		# S-polarized
		imdCff[6, lyrPInd] = /(xfrCffTmpS[1] + xfrCffTmpS[2] * refLS[end], xfrCffTmpS[3] + xfrCffTmpS[4] * refLS[end])
	end
	## Transmission coefficients between cell pairs.
	for srcLyr in srcLyrs
		# Set index for target layer loop. 
		trgLyr = srcLyr + 1
		# Seed expansion construction, note only three coefficients are used in expansion versus 
		# four in nesting.
		# P-polarized.
		trnTmpP .= [1.0 + im * 0.0, 1.0 + im * 0.0, 0.0 + im * 0.0]
		# S-polarized.
		trnTmpS .= [1.0 + im * 0.0, 1.0 + im * 0.0, 0.0 + im * 0.0]
		# Reset target index.
		trgUInd = 1
		# Determine target layers for particular source. 
		trgLyrs = sort(lPairs[2, findall(x -> x == srcLyr, lPairs[1, :])])

		while trgLyr <= trgLyrs[end]

			if trgLyrs[trgUInd] == trgLyr

				for lyrPInd in findall(x -> x == [srcLyr, trgLyr], [lPairs[(2 * ci - 1):(2 * ci)] for ci = 1:size(lPairs)[2]])
					# P-polarized. 
					imdCff[1, lyrPInd] = /(trnTmpP[1] * trnLP[trgLyr - 1], trnTmpP[2] + trnTmpP[3] * refLP[trgLyr - 1])
					# S-polarized. 
					imdCff[2, lyrPInd] = /(trnTmpS[1] * trnLS[trgLyr - 1], trnTmpS[2] + trnTmpS[3] * refLS[trgLyr - 1])
				end
				# Move to next target layer. 
				if length(trgLyrs) > trgUInd
					
					trgUInd += 1
				end
			end
			# Step transmission coefficients. 
			if trgLyr < numLyrs
				# P-polarized
				trnFlmExp!(lVar, enr, wvc, trgLyr, 1, trnLP, trnRP, refLP, refRP, trnTmpP)
				# S-polarized
				trnFlmExp!(lVar, enr, wvc, trgLyr, 1, trnLS, trnRS, refLS, refRS, trnTmpS)
			end
			trgLyr += 1
		end
	end
	# Calculate total field transmission coefficients.
	for lP = 1:size(lPairs)[2]
	
		@views relFlm!(lVar, enr, wvc, lPairs[:,lP], imdCff[:,lP], xfrCff[:,lP])
	end

	return nothing
end
precompile(tfrFlm!, (lyrDsc, Array{Int64,2}, Float64, Float64, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}))
"""
	tfrFunc!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, wvc::Float64, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, scl::Float64, tfrVal::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing

Integrand for heat flux between each layer pair in lPairs, scaled by scl. 
"""
function tfrFunc!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, wvc::Float64, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, scl::Float64, tfrVal::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing
	# Preallocate loop variables. 
	# Permittivity response values. 
	sRsp = 0.0 + im * 0.0
	tRsp = 0.0 + im * 0.0
	# Perpendicular wave vectors.
	swv = 0.0 + im * 0.0
	twv = 0.0 + im * 0.0
	# Layer thickness factors. 
	lTck = 0.0
	sLyrFac = 0.0 
	tLyrFacR = 0.0 
	tLyrFacL = 0.0
	tLyrFacM = 0.0 + im * 0.0
	# Polarization factors. 
	pFacR = 0.0
	pFacL = 0.0
	pFacM = 0.0 + im * 0.0

	for pr = 1 : size(lPairs)[2]
		# Permittivity response. 
		sRsp = lVar.rspPrf[lPairs[1, pr]](enr)
		tRsp = lVar.rspPrf[lPairs[2, pr]](enr)
		# Perpendicular wave vectors.
		swv = wvcPrp(sRsp, lyrTckRel(lVar, lPairs[1, pr], enr), wvc)
		twv = wvcPrp(tRsp, lyrTckRel(lVar, lPairs[2, pr], enr), wvc)
		# Layer thickness factors.
		# Account for possible finite thickness of source layer.
		if (lPairs[1, pr] != 1) && (lPairs[1, pr] != length(lVar.tmpLst))

			lTck = lyrTckRel(lVar, lPairs[1, pr], enr)
			sLyrFac = 1.0 - exp(-4.0 * pi * imag(swv) * lTck)		
		else

			sLyrFac = 1.0
		end
		# Account for possible finite thickness of target layer. 
		# Layer pairs are assumed to be ordered, check made in tfrFlm!.
		if (lPairs[2, pr] != 1) && (lPairs[2, pr] != length(lVar.tmpLst))

			lTck = lyrTckRel(lVar, lPairs[2, pr], enr)
			# Single layer factor comes from shifting factors to avoid overflow.
			tLyrFacR = 1.0 - exp(-4.0 * pi * imag(twv) * lTck)
			tLyrFacL = exp(-4.0 * pi * imag(twv) * lTck) * (1.0 - exp(-4.0 * pi * imag(twv) * lTck))
			tLyrFacM = exp(-4.0 * pi * lTck * (im * real(twv) + imag(twv))) * (1.0 - exp(4.0 * pi * im * real(twv) * lTck))
		else
			
			tLyrFacR = 1.0 
			tLyrFacL = 0.0
			tLyrFacM = 0.0 + im * 0.0
		end
		# Magnitude factor for p-polarized right moving waves.
		pFacR = /(abs(wvc^2 + conj(swv) * twv)^2, abs(tRsp) * abs(sRsp))
		# Magnitude factor for p-polarized left moving waves.
		pFacL = /(abs(wvc^2 - conj(swv) * twv)^2, abs(tRsp) * abs(sRsp))
		# Magnitude factor for p-polarized wave mixing.
		pFacM = /((wvc^2 + conj(swv) * twv) * (wvc^2 - swv * conj(twv)), abs(tRsp) * abs(sRsp))
		## Polarization contributions.
		# P-polarized.
		pPol = /(abs(xfrCff[1, pr])^2 * tLyrFacR * pFacR + abs(xfrCff[2, pr])^2 * tLyrFacL * pFacL, imag(twv)) + /(2.0 * imag(xfrCff[1, pr] * conj(xfrCff[2, pr]) * tLyrFacM * pFacM), real(twv))
		# S-polarized.
		sPol = /(abs(xfrCff[3, pr])^2 * tLyrFacR + abs(xfrCff[4, pr])^2 * tLyrFacL, imag(twv)) + /(2.0 * imag(xfrCff[3, pr] * conj(xfrCff[4, pr]) * tLyrFacM), real(twv))
		# Correct possible NaN underflow.
		if isnan(/(sLyrFac * (pPol + sPol), imag(swv)))

			tfrVal[pr] += 0.0
		else
			# Impose theoretical transfer cut off in case of numerical inaccuracy.
			tfrVal[pr] += scl * wvc * (min(1.0, /(imag(lVar.rspPrf[lPairs[1,pr]](enr)) * imag(lVar.rspPrf[lPairs[2,pr]](enr)) * sLyrFac * pPol, 4.0 * imag(swv) * abs(swv)^2)) + min(1.0, /(imag(lVar.rspPrf[lPairs[1,pr]](enr)) * imag(lVar.tfrFac[lPairs[2,pr]](enr)) * sLyrFac * sPol, 4.0 * imag(swv) * abs(swv)^2)))
		end
	end

	return nothing
end
precompile(tfrFunc!, (lyrDsc, Array{Int64,2}, Float64, Float64, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Float64, Union{Array{Float64,1},SubArray{Float64,1}}))
"""
	evaVarED(dlt::Float64, uvc::Float64)::Float64

Variable transformation for decaying waves. 
"""
function evaVarED(dlt::Float64, uvc::Float64)::Float64

	return tan(/(pi, 2.0) * uvc - dlt)
end
precompile(evaVarED, (Float64, Float64))
"""
	
	tfrFncED!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, uvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, tfrInt::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing

Integrand for evanescent tail of heat flux. 
"""
function tfrFncED!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, uvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, tfrInt::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing
	# Reset value of tfrInt 
	tfrInt .= zeros(Float64, size(lPairs)[2])
	# Calculate field transfer coefficients. 
	tfrFlm!(lVar, lPairs, enr, evaVarED(intDltInf, uvc), bdrArrL, bdrArrR, imdCff, xfrCff)
	# Evanescent integrand for quadrature integration. 
	tfrFunc!(lVar, lPairs, enr, evaVarED(intDltInf, uvc), xfrCff, /(1.0, cos(/(pi, 2.0) * uvc - intDltInf)^2), tfrInt)
	
	return nothing
end
precompile(tfrFncED!, (lyrDsc, Array{Int64,2}, Float64, Float64, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{Float64,1},SubArray{Float64,1}}))
"""
	
	tfrFncED!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, uvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, tfrInt::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing

Threaded integrand for evanescent tail of heat flux. 
"""
function tfrFncEDTH!(lVar::lyrDsc, lPairs::Array{Int64,2}, evlPts::Union{Array{Float64,2},SubArray{Float64,2}}, bdrArrL::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, bdrArrR::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, imdCff::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, xfrCff::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, tfrInt::Union{Array{Float64,2},SubArray{Float64,2}})::Nothing
	
	@threads for evlInd = 1:size(evlPts)[2]
		# Reset value of tfrInt. 
		@views copyto!(tfrInt[:,evlInd], zeros(Float64, size(lPairs)[2]))
		# Calculate field transfer coefficients. 
		@views tfrFlm!(lVar, lPairs, evlPts[1,evlInd], evaVarED(intDltInf, evlPts[2,evlInd]), bdrArrL[:,:,threadid()], bdrArrR[:,:,threadid()], imdCff[:,:,threadid()], xfrCff[:,:,threadid()])
		# Evanescent integrand for quadrature integration. 
		@views tfrFunc!(lVar, lPairs, evlPts[1,evlInd], evaVarED(intDltInf, evlPts[2,evlInd]), xfrCff[:,:,threadid()], /(1.0, cos(/(pi, 2.0) * evlPts[2,evlInd] - intDltInf)^2), tfrInt[:,evlInd])
		# Include heat transfer prefactor J cm^-2
		for lyrPair = 1:size(lPairs)[2]

			tfrInt[lyrPair,evlInd] *= flxPfc(lVar, lPairs[:,lyrPair], evlPts[1,evlInd]) * /(real(lVar.tfrFac[lPairs[1,lyrPair]](evlPts[1,evlInd])) * real(lVar.tfrFac[lPairs[2,lyrPair]](evlPts[1,evlInd])), imag(lVar.rspPrf[lPairs[1,lyrPair]](evlPts[1,evlInd])) * imag(lVar.rspPrf[lPairs[2,lyrPair]](evlPts[1,evlInd])))
		end
	end
	
	return nothing
end
precompile(tfrFncEDTH!, (lyrDsc, Array{Int64,2}, Union{Array{Float64,2},SubArray{Float64,2}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{Float64,2},SubArray{Float64,2}}))
"""
	evaVarEILD(cen::Float64, dlt::Float64, uvc::Float64)::Float64

Double Lorentzian variable transformation for intermediate waves. 
"""
@inline function evaVarEILD(cen::Float64, dlt::Float64, uvc::Float64)::Float64

	return sqrt(/(1.0, sqrt(uvc)) - ^(dlt, 2)) + cen
end
precompile(evaVarEILD, (Float64, Float64, Float64))
"""
	
	tfrFncEILD!(lVar::lyrDsc, lPairs::Array{Int64,2}, cen::Float64, wdt::Float64, dlt::Float64, enr::Float64, uvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, tfrInt::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing

Double Lorentzian scaled integrand for ``intermediate waves''. 
"""
function tfrFncEILD!(lVar::lyrDsc, lPairs::Array{Int64,2}, cen::Float64, wdt::Float64, dlt::Float64, enr::Float64, uvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, tfrInt::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing
	# Reset value of tfrInt 
	tfrInt .= zeros(Float64, size(lPairs)[2])
	# Calculate field transfer coefficients. 
	tfrFlm!(lVar, lPairs, enr, evaVarEILD(cen, dlt, uvc), bdrArrL, bdrArrR, imdCff, xfrCff)
	# Evanescent integrand for quadrature integration. 
	tfrFunc!(lVar, lPairs, enr, evaVarEILD(cen, dlt, uvc), xfrCff, /(1.0, 4.0 * uvc * sqrt(sqrt(uvc) - ^(dlt, 2))), tfrInt)
	# Second portion of integral.
	tfrFlm!(lVar, lPairs, enr, evaVarEILD(cen, dlt, ^(wdt^2 + dlt^2, -2)) - evaVarEILD(cen, dlt, uvc), bdrArrL, bdrArrR, imdCff, xfrCff)
	# Intermediate field integrand for quadrature integration. 
	tfrFunc!(lVar, lPairs, enr, evaVarEILD(cen, dlt, ^(wdt^2 + dlt^2, -2)) - evaVarEILD(cen, dlt, uvc), xfrCff, /(1.0, 4.0 * uvc * sqrt(sqrt(uvc) - ^(dlt, 2))), tfrInt)

	return nothing
end
precompile(tfrFncEILD!, (lyrDsc, Array{Int64,2}, Float64, Float64, Float64, Float64, Float64, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{Float64,1},SubArray{Float64,1}}))
"""
	
	tfrFncPR!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, uvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, tfrInt::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing

Integrand for propagating portion of heat flux. 
"""
function tfrFncPR!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, uvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, xfrCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, tfrInt::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing
	# Reset value of tfrInt. 
	tfrInt .= zeros(Float64, size(lPairs)[2])
	# Calculate field transfer coefficients. 
	tfrFlm!(lVar, lPairs, enr, uvc, bdrArrL, bdrArrR, imdCff, xfrCff)
	# Propagating field integrand for quadrature integration.
	tfrFunc!(lVar, lPairs, enr, uvc, xfrCff, 1.0, tfrInt)

	return nothing
end
precompile(tfrFncPR!, (lyrDsc, Array{Int64,2}, Float64, Float64, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{Float64,1},SubArray{Float64,1}}))
"""
	
	tfrFncPRTH!(lVar::lyrDsc, lPairs::Array{Int64,2}, evlPts::Union{Array{Float64,2},SubArray{Float64,2}}, bdrArrL::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, bdrArrR::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, imdCff::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, xfrCff::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, tfrInt::Union{Array{Float64,2},SubArray{Float64,2}})::Nothing

Threaded integrand for propagating portion of heat flux. 
"""
function tfrFncPRTH!(lVar::lyrDsc, lPairs::Array{Int64,2}, evlPts::Union{Array{Float64,2},SubArray{Float64,2}}, bdrArrL::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, bdrArrR::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, imdCff::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, xfrCff::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, tfrInt::Union{Array{Float64,2},SubArray{Float64,2}})::Nothing
	
	@threads for evlInd = 1:size(evlPts)[2]
		# Reset value of tfrInt. 
		@views copyto!(tfrInt[:,evlInd], zeros(Float64, size(lPairs)[2]))
		# Calculate field transfer coefficients. 
		@views tfrFlm!(lVar, lPairs, evlPts[1,evlInd], evlPts[2,evlInd], bdrArrL[:,:,threadid()], bdrArrR[:,:,threadid()], imdCff[:,:,threadid()], xfrCff[:,:,threadid()])
		# Propagating field integrand for quadrature integration.
		@views tfrFunc!(lVar, lPairs, evlPts[1,evlInd], evlPts[2,evlInd], xfrCff[:,:,threadid()], 1.0, tfrInt[:,evlInd])
		# Include heat transfer prefactor J cm^-2
		for lyrPair = 1:size(lPairs)[2]

			tfrInt[lyrPair,evlInd] *= flxPfc(lVar, lPairs[:,lyrPair], evlPts[1,evlInd]) * /(real(lVar.tfrFac[lPairs[1,lyrPair]](evlPts[1,evlInd])) * real(lVar.tfrFac[lPairs[2,lyrPair]](evlPts[1,evlInd])), imag(lVar.rspPrf[lPairs[1,lyrPair]](evlPts[1,evlInd])) * imag(lVar.rspPrf[lPairs[2,lyrPair]](evlPts[1,evlInd])))
		end
	end

	return nothing
end
precompile(tfrFncPRTH!, (lyrDsc, Array{Int64,2}, Union{Array{Float64,2},SubArray{Float64,2}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{Float64,2},SubArray{Float64,2}}))
"""

	iwLocs(lVar::lyrDsc, enr::Float64)::Array{Float64,2}

Determine location and spacing for intermediate wave transformations. 
Note that output may have NaN entries, and function relies on global constants. 
"""
function iwLocs(lVar::lyrDsc, enr::Float64)::Array{Float64,2}
	
	iwLcs = Float64[]
	iwWdt = Float64[]
	# Find possible points needing more precise treatment.
	for rsp in unique(lVar.rspPrf)

		if /(real(rsp(enr)), imag(rsp(enr))) > iTrnV && imag(rsp(enr)) > 0
			# Locations
			push!(iwLcs, sqrt(real(rsp(enr))))
			# Widths
			push!(iwWdt, imag(rsp(enr)))
		end
	end

	srtInd = sortperm(iwLcs)
	iwLcs .= iwLcs[srtInd]
	iwWdt .= iwWdt[srtInd]
	# Remove locations that are too close together. 
	curr = 1
	
	for ind in 1:(length(iwWdt) - 1)
		
		if ind == 1

			if iwLcs[ind] - iwWdt[ind] < 0.0

				iwWdt[ind] = iwLcs[ind]
			end
		end

		if iwLcs[ind + 1] - iwLcs[curr] < /(iwWdt[curr], 10.0)

			iwLcs[ind + 1] = NaN
			iwWdt[ind + 1] = NaN
		else 

			if /(iwLcs[ind + 1] - iwLcs[curr], 2.0) < iwWdt[curr] 

				iwLcs[curr] = /(iwLcs[ind + 1] - iwLcs[curr], 3.0)
			end

			if /(iwLcs[ind + 1] - iwLcs[curr], 2.0) < iwWdt[ind + 1]

				iwLcs[curr] = /(iwLcs[ind + 1] - iwLcs[curr], 3.0)
			end
			curr = ind + 1
		end
	end
	# Reorder remaining locations. 
	srtInd = sortperm(iwLcs)
	iwLcs .= iwLcs[srtInd]
	iwWdt .= iwWdt[srtInd]
	
	return [iwLcs iwWdt]
end
precompile(iwLocs, (lyrDsc, Float64))
"""
	
	tfrIntEv!(lVar::lyrDsc, lPairs::Array{Int64,2}, enrPts::Union{StepRangeLen{Float64},Array{Float64,1}}, trgPts::Array{Float64,2})::Nothing

Threaded transfer function integrand between all pairs of layers specified in lPairs. 
When multiplied with flxPfc result is heat transfer in J ev^-1 cm^-2. 
"""
function tfrIntEv!(lVar::lyrDsc, lPairs::Array{Int64,2}, enrPts::Union{StepRangeLen{Float64},Array{Float64,1}}, trgPts::Array{Float64,2})::Nothing
	# Number of layer pairs. 
	numLPairs = size(lPairs)[2]
	# Obtain number of active threads. 
	threads = nthreads()
	# Relative wave vector at which to transition to atan variable transformation. 
	wvcDecTrns = 5.0
	# Preallocate working memory for repeated function calls.
	# Boundary field transfer coefficients.
	bdrArrL = Array{ComplexF64,3}(undef, 4, length(lVar.tmpLst) - 1, threads)
	bdrArrR = Array{ComplexF64,3}(undef, 4, length(lVar.tmpLst) - 1, threads)
	# Intermediate layer field coefficients. 
	imdCff = Array{ComplexF64,3}(undef, 10, size(lPairs)[2], threads)
	# Layer pair transfer coefficients. 
	xfrCff = Array{ComplexF64,3}(undef, 4, size(lPairs)[2], threads)
	# Progress bar
	prog = Progress(length(enrPts), dt = 0.5, barglyphs = BarGlyphs("[=> ]"), barlen = 50, color = :green)

	@threads for enrInd = 1 : length(enrPts)
		# Define single variable functions for quadrature integration.
		# Propagating waves
		tfrIntPR!(uvc, intVals) = @views tfrFncPR!(lVar, lPairs, enrPts[enrInd], uvc, bdrArrL[:,:,threadid()], bdrArrR[:,:,threadid()], imdCff[:,:,threadid()], xfrCff[:,:,threadid()], intVals)
		# Rapidly decaying waves, intDltInf sets effective upper bound as tan(pi/2 - inDlt)
		tfrIntED!(uvc, intVals) = @views tfrFncED!(lVar, lPairs, enrPts[enrInd], uvc, bdrArrL[:,:,threadid()], bdrArrR[:,:,threadid()], imdCff[:,:,threadid()], xfrCff[:,:,threadid()], intVals)
		# Treat intermediate waves. 
		# Determine intermediate wave locations.
		iwDsc = iwLocs(lVar, enrPts[enrInd])
		# Account for possible NaN entries.
		nans = findfirst(isnan, view(iwDsc, :, 1))

		if !isnothing(nans)

			numIW = nans - 1
		else

			numIW = length(iwDsc[:,1])
		end
		# Compute integral, filling in portions not treated as intermediate waves with 
		# linear quadrature.
		if iwDsc[1,1] - iwDsc[1,2] > 1.0e-6
		
			trgPts[:,enrInd] .= hquadrature(numLPairs, tfrIntPR!, 0.0, iwDsc[1,1] - iwDsc[1,2], reltol = cubRelTol, error_norm = Cubature.L1)[1]
		else

			trgPts[:,enrInd] .= 0.0
		end

		for ind in 1:numIW

			# Width parameter for double Lorentzian.
			dltIW = /(iwDsc[ind, 2], 10.0)
			# Transformed flux integrand.
			tfrIntEILD!(uvc, intVals) = @views tfrFncEILD!(lVar, lPairs, iwDsc[ind, 1], iwDsc[ind, 2], dltIW, enrPts[enrInd], uvc, bdrArrL[:,:,threadid()], bdrArrR[:,:,threadid()], imdCff[:,:,threadid()], xfrCff[:,:,threadid()], intVals)
			# Perform integrations.
			trgPts[:,enrInd] .+= hquadrature(numLPairs, tfrIntEILD!, ^(^(iwDsc[ind,2], 2) + ^(dltIW, 2), -2), ^(dltIW, -4), reltol = cubRelTol, error_norm = Cubature.L1)[1]

			if ind < numIW
				
				trgPts[:,enrInd] .+= hquadrature(numLPairs, tfrIntPR!, iwDsc[ind, 1] + iwDsc[ind, 2], iwDsc[ind + 1, 1] - iwDsc[ind + 1, 2], reltol = cubRelTol, error_norm = Cubature.L1)[1]
			end	
		end
		# Transition to decaying waves.
		if iwDsc[numIW, 1] + iwDsc[numIW, 2] < wvcDecTrns
		
			trgPts[:,enrInd] .+= hquadrature(numLPairs, tfrIntPR!, iwDsc[numIW, 1] + iwDsc[numIW, 2], wvcDecTrns, reltol = cubRelTol, error_norm = Cubature.L1)[1]
			# Transition point.
			evaDecTrn = /(2.0 * atan(wvcDecTrns), pi)
		else
			# Transition point.
			evaDecTrn = /(2.0 * atan(iwDsc[numIW, 1] + iwDsc[numIW, 2]), pi)
		end
		# Compute decaying wave contribution. 
		trgPts[:,enrInd] .+= hquadrature(numLPairs, tfrIntED!, evaDecTrn, 1.0, reltol = cubRelTol, error_norm = Cubature.L1)[1]
		# Include layer specific prefactor.
		for lyrPair = 1:numLPairs

			trgPts[lyrPair,enrInd] *= flxPfc(lVar, lPairs[:,lyrPair], enrPts[enrInd])
			/(real(lVar.tfrFac[lPairs[1,lyrPair]](enrPts[enrInd])) * real(lVar.tfrFac[lPairs[2,lyrPair]](enrPts[enrInd])),imag(lVar.rspPrf[lPairs[1,lyrPair]](enrPts[enrInd])) * imag(lVar.rspPrf[lPairs[2,lyrPair]](enrPts[enrInd])))
		end
		next!(prog)
	end

	return nothing
end
precompile(tfrIntEv!, (lyrDsc, Array{Int64,2}, Union{StepRangeLen{Float64},Array{Float64,1}}, Array{Float64,2}))
"""
	
	heatTfr!(lVar::lyrDsc, lPairs::Array{Int64,2}, enrRng::Tuple{Float64,Float64}, htPairs::Array{Float64,1})::Nothing

Threaded heat transfer between all pairs of layers specified in lPairs, J cm^-2. 
"""
function heatTfr!(lVar::lyrDsc, lPairs::Array{Int64,2}, enrRng::Tuple{Float64,Float64}, htPairs::Array{Float64,1})::Nothing
	# Number of layer pairs. 
	numLPairs = size(lPairs)[2]
	# Obtain number of active threads. 
	threads = nthreads()
	# Relative wave vector at which to transition to atan variable transformation. 
	wvcDecTrns = 5.0
	## Preallocate working memory for repeated function calls.
	# Boundary field transfer coefficients.
	bdrArrL = Array{ComplexF64,3}(undef, 4, length(lVar.tmpLst) - 1, threads)
	bdrArrR = Array{ComplexF64,3}(undef, 4, length(lVar.tmpLst) - 1, threads)
	# Intermediate layer field coefficients. 
	imdCff = Array{ComplexF64,3}(undef, 10, size(lPairs)[2], threads)
	# Layer pair transfer coefficients. 
	xfrCff = Array{ComplexF64,3}(undef, 4, size(lPairs)[2], threads)
	## Define integrand functions compatible with cubature package. 
	# Propagating waves integrand (no variable redefinition).
	tfrIntPRTH!(evlPts, intVals) = tfrFncPRTH!(lVar, lPairs, evlPts, bdrArrL, bdrArrR, imdCff, xfrCff, intVals)
	# Rapidly decaying waves (atan scaling), intDltInf sets effective upper bound as 
	# tan(pi/2 - inDlt)
	tfrIntEDTH!(evlPts, intVals) = tfrFncEDTH!(lVar, lPairs, evlPts, bdrArrL, bdrArrR, imdCff, xfrCff, intVals)
	## Compute integrals.
	# Up to light line.
	htPairs .= hcubature_v(numLPairs, tfrIntPRTH!, [enrRng[1],0.0], [enrRng[2],1.0], reltol = cubRelTol, error_norm = Cubature.L1)[1]
	# Determine a safer absolute tolerance for further integral contributions.
	cubAbsTol = cubAbsTolCoeff * /(sum(htPairs), length(htPairs))
	# Light line to rapidly decaying transition point.	
	htPairs .+= hcubature_v(numLPairs, tfrIntPRTH!, [enrRng[1],1.0], [enrRng[2],wvcDecTrns], reltol = cubRelTol, abstol = cubAbsTol, error_norm = Cubature.L1)[1]
	# Transform transition point.
	evaDecTrn = /(2.0 * atan(wvcDecTrns), pi)
	# Compute decaying wave contribution. 
	htPairs .+= hcubature_v(numLPairs, tfrIntEDTH!, [enrRng[1],evaDecTrn], [enrRng[2],1.0], reltol = cubRelTol, abstol = cubAbsTol, error_norm = Cubature.L1)[1]
	
	return nothing
end
precompile(heatTfr!, (lyrDsc, Array{Int64,2}, Tuple{Float64,Float64}, Array{Float64,1}))
end 
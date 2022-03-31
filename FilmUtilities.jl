__precompile__()
# See readme for module description. 
module FilmUtilities
using Base.Threads, Cubature
## Load data structures. 
using FilmDataStructures
export plcFnc, subLyrDsc, heatTfr!, relSpont!
## Constants
# Relative tolerance for all cubature calculations.1e-3
const cubRelTol = 1.0e-3
# Sets absolute tolerance for decaying wave cubature relative to (current) calculated propagating part. 1e-2
const cubRatTol = 1.0e-2 
# Relative wavevector at which to transition to tan variable transformation for cubature. 
const wvcDecTrns = 1.0
# Relative range where transformation used for evanescent waves is linear (higher accuracy).
const linRng = 20.0
# Cutoff for unbounded integrals, maximum relative wavevector.
# Should be at least 2 pi lambdaMax / gap, where lambdaMax is the largest relevant value of lambda.
const wvcRelMax = 3.0e5
# Transformation of wavevector cutoff for internal use, edit using wvcRelMax.
const trnWvcMax = /(4.0 * linRng * atan(wvcRelMax), pi) 
# Transformation of wavevector transition point for internal use, edit using wvcDecTrns.
const trnWvcMin = /(4.0 * linRng * atan(wvcDecTrns), pi) 
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
Analytic extension of the geometric series for a given initial term and geometric factor.
"""
@inline function geoSrs(intTrm::ComplexF64, geoFct::ComplexF64)::ComplexF64

	return /(intTrm, 1.0 - geoFct)
end
precompile(geoSrs, (ComplexF64, ComplexF64))
"""
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
Normalized perpendicular wave, incorporating an effective normalization.
"""
@inline function wvcPrp(rsp::ComplexF64, wvc::Float64)::ComplexF64

	return zSqrt(rsp - ^(wvc, 2.0)) 
end
precompile(wvcPrp, (ComplexF64, Float64))
"""

	plcFnc(enr::Float64, tmp::Float64)::Float64

Planck function for energy in electron volts and temperature in Kelvin.
"""
@inline function plcFnc(enr::Float64, tmp::Float64)::Float64

	return /(enr, exp(/(enr, blzK * tmp)) - 1.0)
end
precompile(plcFnc, (Float64, Float64))
"""
Power prefactor, W eV^-1 cm^-2, for fluctuating electric field intensity integrand. 
"""
@inline function flxPfc(lVar::lyrDsc, lPair::Union{Array{Int64,1},SubArray{Int64,1}}, enr::Float64)::Float64

	return /(evJ * enr^2, muEv^2 * hBEv * ^(10.0, -8)) * (plcFnc(enr, lVar.tmpLst[lPair[1]]) - plcFnc(enr, lVar.tmpLst[lPair[2]]))
end
precompile(flxPfc, (lyrDsc, Union{Array{Int64,1},SubArray{Int64,1}}, Float64))
"""
Calculation of phase accumulation under propagation, distance assumed to be in units of wavelength.
"""
@inline function prpLyr(dst::Float64, rsp::ComplexF64, wvc::Float64)::ComplexF64

	return exp(2.0 * pi * im * wvcPrp(rsp, wvc) * dst)
end
precompile(prpLyr, (Float64, ComplexF64, Float64))
"""
Layer thickness in microns. 
"""
@inline function lyrTck(lVar::lyrDsc, lNum::Int)::Float64

	if lNum == 1 || lNum == length(lVar.tmpLst)

		return 0.0
	else

		return /(lVar.bdrLoc[lNum] - lVar.bdrLoc[lNum - 1], 1.0)
	end
end
precompile(lyrTck, (lyrDsc, Int))
"""
Compute layer thickness in units of relative wavelength.
"""
@inline function lyrTckRel(lVar::lyrDsc, lNum::Int, enr::Float64)::Float64

	return /(enr * lyrTck(lVar, lNum), muEv)
end
precompile(lyrTckRel, (lyrDsc, Int, Float64))
"""
Compute total structure thickness in units of relative wavelength.
"""
@inline function strTck(lVar::lyrDsc, enr::Float64)::Float64

	return /(enr * (lVar.bdrLoc[end] - lVar.bdrLoc[1]), muEv)
end
precompile(strTck, (lyrDsc, Float64))
"""
Isotopic electric reflection and transmission coefficients for two half spaces. 
First entries in the permittivity (rps) and wave vector (kPrp) tuples are consider to be 
incident. 
The first two entries returned are for 'p' polarized waves, the next two for 's' 
polarized waves. 
Odd numbers are transmission, even modes are reflection.
"""
@inline function xfrBdr!(kPrp::NTuple{2,ComplexF64}, rsp::NTuple{2,ComplexF64}, bdrArr::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing
	# Fresnel coefficients.
	copyto!(bdrArr, [/(2.0 * zSqrt(rsp[1] * rsp[2]) * kPrp[1], rsp[1] * kPrp[2] + rsp[2] * kPrp[1]), /(rsp[1] * kPrp[2] - rsp[2] * kPrp[1], rsp[1] * kPrp[2] + rsp[2] * kPrp[1]), /(2.0 * kPrp[1], kPrp[2] + kPrp[1]), /(kPrp[2] - kPrp[1], kPrp[2] + kPrp[1])])

	 return nothing
end
precompile(xfrBdr!, (NTuple{2,ComplexF64}, NTuple{2,ComplexF64}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}))
"""
Fills bdrArr with individual boundary transfer coefficients, i.e. as though every boundary were 
between half-spaces. 
If mode == 1 smaller indices are treated as initial. 
If mode == 2 larger indices are treated as initial.
Convention for a given column of bdrArr as described for xfrBdr. 
"""
@inline function bdrFill!(lVar::lyrDsc, enr::Float64, wvc::Float64, mode::Int, bdrArr::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}})::Nothing
	# Smaller indices are treated as initial. 
	if mode == 1

		for bdr = 1 : length(lVar.bdrLoc)

			@views xfrBdr!((wvcPrp(lVar.rspPrf[bdr](enr), wvc), wvcPrp(lVar.rspPrf[bdr + 1](enr), wvc)), (lVar.rspPrf[bdr](enr), lVar.rspPrf[bdr + 1](enr)), bdrArr[:,bdr])
		end
	# Larger indices are treated as initial. 
	elseif mode == 2

		for bdr = 1 : length(lVar.bdrLoc)
			
			@views xfrBdr!((wvcPrp(lVar.rspPrf[bdr + 1](enr), wvc), wvcPrp(lVar.rspPrf[bdr](enr), wvc)), (lVar.rspPrf[bdr + 1](enr), lVar.rspPrf[bdr](enr)), bdrArr[:,bdr])
		end
	else

		error("Unrecognized fill mode.") 
		return nothing	
	end

	return nothing 
end
precompile(bdrFill!, (lyrDsc, Float64, Float64, Int, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}))
"""
Updates transmission coefficients under expansion of nested reflection coefficient. 
See associated notes for additional details. 
If mode == 1, expansion is ``to the right'', i.e. expanding to growing layers numbers. 
If mode == 2, expansion is ``to the left'', i.e. expanding to shrinking layers numbers.
First entry of trnArr should be numerator coefficient. 
The next two entries are the reflection and additional coefficients for the denominator. 
See notes for additional details. 
"""
@inline function trnFlmExp!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, trnArrHL::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, trnArrHR::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, refArrHL::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, refArrHR::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, trnArr::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing

	if mode == 1
	
		copyto!(trnArr, [trnArrHL[lNum - 1] * prpLyr(lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * trnArr[1], trnArr[2] + refArrHL[lNum - 1] * trnArr[3], prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * (trnArr[3] - refArrHR[lNum - 1] * trnArr[2])])

	elseif mode == 2
		
		copyto!(trnArr, [trnArrHR[lNum] * prpLyr(lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * trnArr[1], trnArr[2] + refArrHR[lNum] * trnArr[3], prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum ](enr), wvc) * (trnArr[3] - refArrHL[lNum] * trnArr[2])])
	else 

		error("Unrecognized increment mode.")	
		return nothing
	end

	return nothing
end
precompile(trnFlmExp!, (lyrDsc, Float64, Float64, Int, Int, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}))
"""
Updates transfer coefficients by nesting current values, i.e. inclusion of an additional layer. 
See associated notes for additional details. 
If mode == 1, nesting is ``to the right'', i.e. adding growing layer numbers. 
If mode == 2, nesting is ``to the left'', i.e. adding shrinking layer numbers.
xfrArr should hold the expansion coefficients, to the current level, of the numerator as  
its first two entries, additional term followed by reflection term. 
The next two entries follow the same convention, but for the denominator. 
"""
@inline function xfrFlmNst!(lVar::lyrDsc, enr::Float64, wvc::Float64, lNum::Int, mode::Int, refArrHL::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, refArrHR::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, xfrArr::Union{Array{ComplexF64,1},SubArray{ComplexF64,1}})::Nothing

	if mode == 1
	
		copyto!(xfrArr, [refArrHR[lNum] * xfrArr[3] + prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * xfrArr[1], refArrHR[lNum] * xfrArr[4] + prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * xfrArr[2], xfrArr[3] - prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * refArrHL[lNum] * xfrArr[1], xfrArr[4] - prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * refArrHL[lNum] * xfrArr[2]])
	elseif mode == 2
		
		copyto!(xfrArr, [refArrHL[lNum - 1] * xfrArr[3] + prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * xfrArr[1], refArrHL[lNum - 1] * xfrArr[4] + prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * xfrArr[2], xfrArr[3] - prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * refArrHR[lNum - 1] * xfrArr[1], xfrArr[4] - prpLyr(2.0 * lyrTckRel(lVar, lNum, enr), lVar.rspPrf[lNum](enr), wvc) * refArrHR[lNum - 1] * xfrArr[2]])
	else 

		error("Unrecognized increment mode.")	
		return nothing
	end

	return nothing
end
precompile(xfrFlmNst!, (lyrDsc, Float64, Float64, Int, Int, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}, Union{Array{ComplexF64,1},SubArray{ComplexF64,1}}))
"""
Calculates field transfer coefficients for a list of ordered pairs of layers. 
bdrArrL will contain interface transfer coefficients with smaller layer numbers treated as 
initial. 
bdrArrR will contain interface transfer coefficients with larger layer numbers treated as 
initial. 
"""
function tfrFlm!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, wvc::Float64, bdrArrL::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, bdrArrR::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}})::Nothing
	# Calculate infinite interface coefficients.
	# Left incident.
	bdrFill!(lVar, enr, wvc, 1, bdrArrL)
	# Right incident.
	bdrFill!(lVar, enr, wvc, 2, bdrArrR)
	# Double check that pairs in list are correctly ordered. 
	for pr = 1:size(lPairs)[2]

		if lPairs[2, pr] < lPairs[1, pr]
			
			error("Current program convention requires ordered layers. Transfer disregarding flux prefactor for heat transfer---mode = 1---is symmetric under pair reversal.")
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
		ind += 1
	end
	# Save right incident reflection coefficients for largest indexed source cells (if coinciding with target).
	if srcLyrs[end] == trgLyrs[end]
		
		for lyrPInd in findall(x -> x == srcLyrs[end], lPairs[1, :])
			# P-polarized.
			imdCff[3, lyrPInd] = /(xfrCffTmpP[1] + xfrCffTmpP[2] * refRP[1], xfrCffTmpP[3] + xfrCffTmpP[4] * refRP[1])
			# S-polarized.
			imdCff[4, lyrPInd] = /(xfrCffTmpS[1] + xfrCffTmpS[2] * refRS[1], xfrCffTmpS[3] + xfrCffTmpS[4] * refRS[1])
		end
	end
	# Save right incident reflection coefficients for largest indexed target cells.
	for lyrPInd in findall(x -> x == trgLyrs[end], lPairs[2, :])
		# P-polarized.
		imdCff[7, lyrPInd] = /(xfrCffTmpP[1] + xfrCffTmpP[2] * refRP[1], xfrCffTmpP[3] + xfrCffTmpP[4] * refRP[1])
		# S-polarized.
		imdCff[8, lyrPInd] = /(xfrCffTmpS[1] + xfrCffTmpS[2] * refRS[1], xfrCffTmpS[3] + xfrCffTmpS[4] * refRS[1])
	end
	## Left incident reflection coefficients.
	# Safeguard against fictitious reflections from infinite boundary for targets.
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
	return nothing
end
precompile(tfrFlm!, (lyrDsc, Array{Int64,2}, Float64, Float64, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}))
"""
Integrand for heat flux between each layer pair in lPairs, scaled by scl. 
imdCff[:, lyr] = [trnCf; refLeftCellLeft; refRightCellLeft; refLeftCellRight; 
refRightCellRight], with p-pol followed by s-pol in each. 
"""
function tfrFunc!(lVar::lyrDsc, lPairs::Array{Int64,2}, enr::Float64, wvc::Float64, imdCff::Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, scl::Float64, tfrVal::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing
	# Preallocate loop variables. 
	# Permittivity response values. 
	sRsp = 0.0 + im * 0.0
	tRsp = 0.0 + im * 0.0
	# Perpendicular wave vectors.
	swv = 0.0 + im * 0.0
	twv = 0.0 + im * 0.0
	# Layer thickness factors. 
	slTck = 0.0
	tlTck = 0.0
	sLyrFacR = 0.0 
	sLyrFacL = 0.0
	sLyrFacM = 0.0 + im * 0.0
	tLyrFacR = 0.0 
	tLyrFacL = 0.0
	tLyrFacM = 0.0 + im * 0.0
	# Polarization factors. 
	pFacW = 0.0
	pFacM = 0.0 
	# Layer propagation variables.
	lLyrPhz = 0.0 + im * 0.0
	rLyrPhz = 0.0 + im * 0.0
	# Wave dressing factors.
	pPolSrc = 0.0 
	sPolSrc = 0.0
	pPolTrg = 0.0
	sPolTrg = 0.0
	# Geometric dressing factors.
	pGeo = 0.0
	sGeo = 0.0

	for pr = 1 : size(lPairs)[2]
		# Permittivity response. 
		sRsp = lVar.rspPrf[lPairs[1, pr]](enr)
		tRsp = lVar.rspPrf[lPairs[2, pr]](enr)
		# Perpendicular wave vectors.
		swv = wvcPrp(sRsp, wvc)
		twv = wvcPrp(tRsp, wvc)
		
		if lPairs[1, pr] == lPairs[2, pr] 

			error("Can not compute self layer transfer---answer is not finite.")
			return nothing
		end	
		# Account for possible finite thickness of source layer.
		if (lPairs[1, pr] != 1) && (lPairs[1, pr] != length(lVar.tmpLst))

			# Layer thickness factors.
			slTck = lyrTckRel(lVar, lPairs[1, pr], enr)
			sLyrFacR = 1.0 - exp(-4.0 * pi * abs(imag(swv)) * slTck)
			sLyrFacL = exp(-4.0 * pi * abs(imag(swv)) * slTck) * (1.0 - exp(-4.0 * pi * abs(imag(swv)) * slTck))
			sLyrFacM = exp(-4.0 * pi * abs(imag(swv)) * slTck) * (1.0 - exp(4.0 * pi * im * real(swv) * slTck))
			# Layer phase propagation
			lLyrPhz = prpLyr(slTck, lVar.rspPrf[lPairs[1, pr]](enr), wvc)
		else

			sLyrFacR = 1.0 
			sLyrFacL = 0.0
			sLyrFacM = 0.0 + im * 0.0
			lLyrPhz = 0.0
		end
		# Account for possible finite thickness of target layer. 
		# Layer pairs are assumed to be ordered, check made in tfrFlm!.
		if (lPairs[2, pr] != 1) && (lPairs[2, pr] != length(lVar.tmpLst))

			# Layer thickness factors.
			tlTck = lyrTckRel(lVar, lPairs[2, pr], enr)
			# Single layer factor comes from shifting factors to avoid overflow.
			tLyrFacR = 1.0 - exp(-4.0 * pi * abs(imag(twv)) * tlTck)
			tLyrFacL = exp(-4.0 * pi * abs(imag(twv)) * tlTck) * (1.0 - exp(-4.0 * pi * abs(imag(twv)) * tlTck))
			tLyrFacM = exp(-4.0 * pi * abs(imag(twv)) * tlTck) * (1.0 - exp(4.0 * pi * im * real(twv) * tlTck))
			# Layer phase propagation
			rLyrPhz = prpLyr(tlTck, lVar.rspPrf[lPairs[2, pr]](enr), wvc)
		else
			
			tLyrFacR = 1.0 
			tLyrFacL = 0.0
			tLyrFacM = 0.0 + im * 0.0
			rLyrPhz = 0.0

		end
		# Magnitude factor for p-polarized for waves moving in same direction. 
		pFacW = /((^(wvc, 2.0) + ^(abs(swv), 2.0)) * (^(wvc, 2.0) + ^(abs(twv), 2.0)), abs(tRsp) * abs(sRsp))
		# Magnitude factor for p-polarized wave mixing.
		pFacM = /((^(wvc, 2.0) - ^(abs(swv), 2.0)) * (^(wvc, 2.0) - ^(abs(twv), 2.0)), abs(tRsp) * abs(sRsp))
		## Wave mixing (left-right) contributions in source layer.
		# P-polarized.
		pPolSrc = /(sLyrFacR + abs(imdCff[3, pr])^2 * sLyrFacL, imag(swv)) + /(imag(2.0 * conj(imdCff[3, pr]) * sLyrFacM), real(swv)) * /(pFacM, pFacW)
		# S-polarized.
		sPolSrc = /(sLyrFacR + abs(imdCff[4, pr])^2 * sLyrFacL, imag(swv)) + /(imag(2.0 * conj(imdCff[4, pr]) *  sLyrFacM), real(swv)) 
		## Wave mixing (left-right) contributions in target layer.
		# P-polarized.
		pPolTrg = /(tLyrFacR + abs(imdCff[9, pr])^2 * tLyrFacL, imag(twv)) + /(imag(2.0 * conj(imdCff[9, pr]) * tLyrFacM), real(twv)) * /(pFacM, pFacW)
		# S-polarized.
		sPolTrg = /(tLyrFacR + abs(imdCff[10, pr])^2 * tLyrFacL, imag(twv)) + /(imag(2.0 * conj(imdCff[10, pr]) * tLyrFacM), real(twv))
		## Layer-to-layer transfer coefficients.
		# P-polarized.
		pPolTrf = ^(abs(/(imdCff[1, pr], 2.0 * swv)), 2.0)
		# S-polarized.
		sPolTrf = ^(abs(/(imdCff[2, pr], 2.0 * swv)), 2.0)
		## Geometric dressing factors, distinct from field mixing in the source and target layers.
		# P-polarized.
		pGeo = ^(abs(geoSrs(1.0 + im * 0.0, imdCff[3, pr] * lLyrPhz * imdCff[5, pr] * lLyrPhz) * geoSrs(1.0 + im * 0.0, imdCff[7, pr] * rLyrPhz * imdCff[9, pr] * rLyrPhz)), 2.0)
		# S-polarized.
		sGeo = ^(abs(geoSrs(1.0 + im * 0.0, imdCff[4, pr] * lLyrPhz * imdCff[6, pr] * lLyrPhz) * geoSrs(1.0 + im * 0.0, imdCff[8, pr] * rLyrPhz * imdCff[10, pr] * rLyrPhz)), 2.0)
		# Correct possible NaN underflow.
		if isnan(/(abs(imdCff[1, pr])^2 * pPolSrc * pPolTrg * pGeo * pFacW + abs(imdCff[2, pr])^2 * sPolSrc * sPolTrg * sGeo, imag(swv)))
			
		else
			# Impose theoretical transfer coefficient cut off in case of numerical inaccuracy.
			tfrVal[pr] += scl * wvc * (min(1.0, imag(lVar.rspPrf[lPairs[1,pr]](enr)) * lVar.tfrFac[lPairs[2,pr]](enr) * pPolTrf * pPolSrc * pPolTrg * pGeo * pFacW) + min(1.0, imag(lVar.rspPrf[lPairs[1,pr]](enr)) * lVar.tfrFac[lPairs[2,pr]](enr) * sPolTrf * sPolSrc * sPolTrg * sGeo))
		end
	end

	return nothing
end
precompile(tfrFunc!, (lyrDsc, Array{Int64,2}, Float64, Float64, Union{Array{ComplexF64,2},SubArray{ComplexF64,2}}, Float64, Union{Array{Float64,1},SubArray{Float64,1}}))
"""
Variable transformation for decaying waves. 
"""
@inline function evaVarED(uvc::Float64)::Float64

	return tan(/(pi, 4.0 * linRng) * uvc )
end
precompile(evaVarED, (Float64, ))
"""
Form scaling for decaying wave transformation defined by evaVarED.
"""
@inline function evaFrmED(uvc::Float64)::Float64

	return /(pi, 4.0 * linRng) * ^(cos(/(pi * uvc, 4.0 * linRng)), -2.0)
end
precompile(evaFrmED, (Float64, ))
"""
Threaded integrand for evanescent tail of electromagnetic transfer function.
mode == 1 is bare transfer function, mode == 2 includes heat flux prefactor. 
"""
function tfrFncEDTH!(mode::Int, lVar::lyrDsc, lPairs::Array{Int64,2}, evlPts::Union{Array{Float64,2},SubArray{Float64,2}}, bdrArrL::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, bdrArrR::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, imdCff::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, tfrInt::Union{Array{Float64,2},SubArray{Float64,2}})::Nothing
	
	@threads for evlInd = 1:size(evlPts)[2]
		# Reset value of tfrInt. 
		@views copyto!(tfrInt[:,evlInd], zeros(Float64, size(lPairs)[2]))
		# Calculate field transfer coefficients. 
		@views tfrFlm!(lVar, lPairs, evlPts[1,evlInd], evaVarED(evlPts[2,evlInd]), bdrArrL[:,:,threadid()], bdrArrR[:,:,threadid()], imdCff[:,:,threadid()])
		# Evanescent integrand for quadrature integration. 
		@views tfrFunc!(lVar, lPairs, evlPts[1,evlInd], evaVarED(evlPts[2,evlInd]), imdCff[:,:,threadid()], evaFrmED(evlPts[2,evlInd]), tfrInt[:,evlInd])
		# Include heat transfer prefactor J cm^-2
		for lyrPair = 1:size(lPairs)[2]

			if mode == 1
			
			elseif 	mode == 2
				
				tfrInt[lyrPair,evlInd] *= flxPfc(lVar, lPairs[:,lyrPair], evlPts[1,evlInd]) 
			else

				error("Unrecognized operation mode.")	
				return nothing
			end
		end
	end
	
	return nothing
end
precompile(tfrFncEDTH!, (Int, lyrDsc, Array{Int64,2}, Union{Array{Float64,2},SubArray{Float64,2}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{Float64,2},SubArray{Float64,2}}))
"""
	
	tfrFncPRTH!(mode::Int, lVar::lyrDsc, lPairs::Array{Int64,2},
	evlPts::Union{Array{Float64,2},SubArray{Float64,2}}, 
	bdrArrL::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, 
	bdrArrR::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, 
	imdCff::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, 
	tfrInt::Union{Array{Float64,2},SubArray{Float64,2}})::Nothing

Threaded integrand for propagating portion for electromagnetic transfer 
function. mode == 1 is bare transfer function, mode == 2 includes heat 
flux prefactor. 
"""
function tfrFncPRTH!(mode::Int, lVar::lyrDsc, lPairs::Array{Int64,2}, evlPts::Union{Array{Float64,2},SubArray{Float64,2}}, bdrArrL::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, bdrArrR::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, imdCff::Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, tfrInt::Union{Array{Float64,2},SubArray{Float64,2}})::Nothing
	
	@threads for evlInd = 1:size(evlPts)[2]
		# Reset value of tfrInt. 
		@views copyto!(tfrInt[:,evlInd], zeros(Float64, size(lPairs)[2]))
		# Calculate field transfer coefficients. 
		@views tfrFlm!(lVar, lPairs, evlPts[1,evlInd], evlPts[2,evlInd], bdrArrL[:,:,threadid()], bdrArrR[:,:,threadid()], imdCff[:,:,threadid()])
		# Propagating field integrand for quadrature integration.
		@views tfrFunc!(lVar, lPairs, evlPts[1,evlInd], evlPts[2,evlInd], imdCff[:,:,threadid()], 1.0, tfrInt[:,evlInd])
		# Include heat transfer prefactor W cm^-2 eV^-1
		for lyrPair = 1:size(lPairs)[2]

			if mode == 1
			
			elseif 	mode == 2
				
				tfrInt[lyrPair,evlInd] *= flxPfc(lVar, lPairs[:,lyrPair], evlPts[1,evlInd]) 
			else

				error("Unrecognized operation mode.")
				return nothing
			end
		end
	end

	return nothing
end
precompile(tfrFncPRTH!, (Int, lyrDsc, Array{Int64,2}, Union{Array{Float64,2},SubArray{Float64,2}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{ComplexF64,3},SubArray{ComplexF64,3}}, Union{Array{Float64,2},SubArray{Float64,2}}))
"""
	
	heatTfr!(lVar::lyrDsc, lPairs::Array{Int64,2}, 
	enrRng::Tuple{Float64,Float64}, htPairs::Array{Float64,1})::Nothing

Threaded calculation of heat transfer between all pairs of layers 
specified in lPairs, reported in units of W cm^-2. 
"""
function heatTfr!(lVar::lyrDsc, lPairs::Array{Int64,2}, enrRng::Tuple{Float64,Float64}, htPairs::Array{Float64,1})::Nothing
	# Number of layer pairs. 
	numLPairs = size(lPairs)[2]
	# Obtain number of active threads. 
	threads = nthreads()
	## Preallocate working memory for repeated function calls.
	# Boundary field transfer coefficients.
	bdrArrL = Array{ComplexF64,3}(undef, 4, length(lVar.tmpLst) - 1, threads)
	bdrArrR = Array{ComplexF64,3}(undef, 4, length(lVar.tmpLst) - 1, threads)
	# Intermediate layer field coefficients. 
	imdCff = Array{ComplexF64,3}(undef, 10, size(lPairs)[2], threads)
	## Define integrand functions compatible with cubature package. 
	# Propagating waves integrand (no variable redefinition).
	tfrIntPRTH!(evlPts, intVals) = tfrFncPRTH!(2, lVar, lPairs, evlPts, bdrArrL, bdrArrR, imdCff, intVals)
	# Rapidly decaying waves (tan scaling)
	tfrIntEDTH!(evlPts, intVals) = tfrFncEDTH!(2, lVar, lPairs, evlPts, bdrArrL, bdrArrR, imdCff, intVals)
	## Compute integrals.
	# Up to light line.
	htPairs .= hcubature_v(numLPairs, tfrIntPRTH!, [enrRng[1], 0.0], [enrRng[2], 1.0], reltol = cubRelTol, error_norm = Cubature.L1)[1]
	# Determine a safer absolute tolerance for further integral contributions.
	cubAbsTol = cubRatTol * /(sum(htPairs), length(htPairs))
	# Light line to rapidly decaying transition point.	
	htPairs .+= hcubature_v(numLPairs, tfrIntPRTH!, [enrRng[1], 1.0], [enrRng[2], wvcDecTrns], reltol = cubRelTol, abstol = cubAbsTol, error_norm = Cubature.L1)[1]
	# Compute decaying wave contribution, range transforms already taken care of. 
	htPairs .+= hcubature_v(numLPairs, tfrIntEDTH!, [enrRng[1], trnWvcMin], [enrRng[2], trnWvcMax], reltol = cubRelTol, abstol = cubAbsTol, error_norm = Cubature.L1)[1]
	return nothing
end
precompile(heatTfr!, (lyrDsc, Array{Int64,2}, Tuple{Float64,Float64}, Array{Float64,1}))
"""

	subLyrDsc(lVar::lyrDsc, lyrNum::Int, lyrDiv::Int)::lyrDsc

Subdivides layer of a given layer description, creating new description.
"""
@inline function subLyrDsc(lVar::lyrDsc, lyrNum::Int, lyrDiv::Int)::lyrDsc

	return lyrDsc([lVar.bdrLoc[1:(lyrNum - 1)]; LinRange(lVar.bdrLoc[lyrNum - 1], lVar.bdrLoc[lyrNum], lyrDiv + 1)[2:(end - 1)]; lVar.bdrLoc[lyrNum:end]], [lVar.tmpLst[1:(lyrNum - 1)]; repeat([lVar.tmpLst[lyrNum]], lyrDiv); lVar.tmpLst[(lyrNum + 1):end]], [lVar.rspPrf[1:(lyrNum - 1)]; repeat([lVar.rspPrf[lyrNum]], lyrDiv); lVar.rspPrf[(lyrNum + 1):end]], [lVar.tfrFac[1:(lyrNum - 1)]; repeat([lVar.tfrFac[lyrNum]], lyrDiv); lVar.tfrFac[(lyrNum + 1):end]])
end
precompile(subLyrDsc, (lyrDsc, Int, Int))
"""

	relSpont!(lVar::lyrDsc, lyrNum::Int, lyrDiv::Int, 
	enrPts::Union{StepRangeLen{Float64},Array{Float64,1}}, 
	relSpont::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing

Approximate threaded calculation of relative spontaneous emission rate 
at any layer in slab structure. To use this function, all / any vacuum 
layers must be given a small imaginary component. Cut-off is set by 
removing self-layer computation. 
"""
function relSpont!(lVar::lyrDsc, lyrNum::Int, enrRng::Tuple{Float64,Float64}, relSpont::Union{Array{Float64,1},SubArray{Float64,1}})::Nothing
	
	if (lyrNum == 1) || (lyrNum == length(lVar.tmpLst))
		
		error("Can not compute relative spontaneous emission rate of infinite layer---division is not sensible. To calculate the relative spontaneous emission rate at this position, create a sub-divided layer description.")
		return nothing
	end	
	# Obtain number of active threads. 
	threads = nthreads()
 	## Create layer pair variable.	
 	numLyrs = length(lVar.tmpLst)
 	numLyrPairs = numLyrs - 1
	lPairs = Array{Int64,2}(undef, 2, numLyrPairs)
	# Fill layer pairs.
	prInd = 1

	for evalLyr = 1 : numLyrs

		if lyrNum < evalLyr 

			lPairsD[1,prInd] = spontLyr
			lPairsD[2,prInd] = evalLyr
			prInd += 1

		elseif lyrNum > evalLyr

			lPairsD[1,prInd] = evalLyr
			lPairsD[2,prInd] = spontLyr
			prInd += 1
		else
		end
	end
	## Preallocate working memory for repeated function calls.
	# Boundary field transfer coefficients.
	bdrArrL = Array{ComplexF64,3}(undef, 4, numLyrPairs, threads)
	bdrArrR = Array{ComplexF64,3}(undef, 4, numLyrPairs, threads)
	# Intermediate layer field coefficients. 
	imdCff = Array{ComplexF64,3}(undef, 10, numLyrPairs, threads)
	# Temporary container for transfer values. 
	trfPairs = Array{Float64,1}(undef, numLyrPairs)
	## Define integrand functions compatible with cubature package. 
	# Propagating waves integrand (no variable redefinition).
	tfrIntPRTH!(evlPts, intVals) = tfrFncPRTH!(1, lVar, lPairs, evlPts, bdrArrL, bdrArrR, imdCff, intVals)
	# Rapidly decaying waves (tan scaling)
	tfrIntEDTH!(evlPts, intVals) = tfrFncEDTH!(1, lVar, lPairs, evlPts, bdrArrL, bdrArrR, imdCff, intVals)
	## Compute pair integrals.
	# Up to light line.
	trfPairs .= hcubature_v(numLyrPairs, tfrIntPRTH!, [enrRng[1], 0.0], [enrRng[2], 1.0], reltol = cubRelTol, error_norm = Cubature.L1)[1]
	# Determine a safer absolute tolerance for further integral contributions.
	cubAbsTol = cubRatTol * /(sum(trfPairs), length(trfPairs))
	# Light line to rapidly decaying transition point.	
	trfPairs .+= hcubature_v(numLyrPairs, tfrIntPRTH!, [enrRng[1], 1.0], [enrRng[2], wvcDecTrns], reltol = cubRelTol, abstol = cubAbsTol, error_norm = Cubature.L1)[1]
	# Compute decaying wave contribution, range transforms already taken care of. 
	trfPairs .+= hcubature_v(numLyrPairs, tfrIntEDTH!, [enrRng[1], trnWvcMin], [enrRng[2], trnWvcMax], reltol = cubRelTol, abstol = cubAbsTol, error_norm = Cubature.L1)[1]
	# Sum results.
	for subLyr = 1 : lyrDiv

		relSpont[1] = 0.0
	end

	for subLyr = 1 : lyrDiv, prtn = 1 : (numLyrD - 1)

		if isnan(trfPairs[(subLyr - 1) * (numLyrD - 1) + prtn])

		else
			
			relSpont[1] += trfPairs[(subLyr - 1) * (numLyrD - 1) + prtn]
		end
	end
	## Normalize spontaneous emission results.
	return nothing
end
precompile(relSpont!, (lyrDsc, Int, Int, Tuple{Float64,Float64}, Union{Array{Float64,1},SubArray{Float64,1}}))
end 
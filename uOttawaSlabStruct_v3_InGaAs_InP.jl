# Construction of PV structure for University of Ottawa collaboration.
# See slides for schematic, slab structure is Si,Gap, Si, Gap, InP, InGaAs, InGaAs, InP, gold.
# All temperatures are in Kelvin, and all lengths are in microns.
# Absorption functions are current setup to calculate the number of photons absorbed as interband transitions.
# - v3: doping for FSF layer, counts above interband photon absorption numbers 
#Input: Heat_or_Photon , 0 == Heat , 1 == Photon
####


const e = 1.602*^(10,-19)
const eV = 1.60218*^(10,-19)
const hbar = 1.05457*^(10,-34) #m2kg/s
const hbEV = 6.5821*^(10,-16) #eV s
function uOttawaSlabs_v3_InGaAs_InP(tEmit::Float64, tBck::Float64, divProtCell::Int, divNCell::Int, divPCell::Int, divSubCell::Int, firstGap::Float64, thckRad::Float64, distGap::Float64, thckProt::Float64, thckInAs::Float64, thckInAsSbP::Float64, thckSub::Float64, Si_ndoping::Float64, prot_ndoping::Float64, ndoping_InAs::Float64, pdoping_InAs::Float64, substrate_doping::Float64, quat_x::Float64, prot_quat_x::Float64, mboxProt::Float64, mboxN::Float64, mboxP::Float64, mboxSub::Float64,xMinProt::Float64, xMinN::Float64, xMinP::Float64, xMinSub::Float64,Heat_or_Photon::Int)
### Settings
	# Total number of layers.
	if firstGap == 0.0
		numExtraLayers = 3
	else
		numExtraLayers = 4
	end
	numLayers = divProtCell + divNCell + divPCell + divSubCell + numExtraLayers
	## Temperature list of the layers.
	tmpLst = fill(tBck, numLayers)
	# Set temperature of the emitter.
	tmpLst[numExtraLayers-2] = tEmit
	## Boundary locations.
	bdrLoc = Vector{Float64}(undef, numLayers - 1) 
	bdrLoc[1] = 0.0
	if firstGap == 0.0
		bdrLoc[2] = distGap
	else
		bdrLoc[2] = thckRad  
		bdrLoc[3] = distGap + thckRad
	end
	#fill boundaries with protection layer, using the graded layer thcikness scheme
	for ind = 1 : divProtCell
		bdrLoc[numExtraLayers - 1 + ind] = bdrLoc[numExtraLayers - 2  + ind] + thckProt/divProtCell
	end
	#=
	for ind = 0 : divProtCell -1
		bdrLoc[4 + ind + 1] = bdrLoc[4 + ind] + xMinProt*(mboxProt^ind)
	end
	=#
	# Fill boundaries for N-type part of the cell.
	#=
	for ind = 0 : divNCell - 1
		bdrLoc[4 + divProtCell + ind + 1] = bdrLoc[4 + divProtCell + ind] + xMinN*(mboxN^ind)
	end
	=#
	for ind = 1 : divNCell

		bdrLoc[numExtraLayers - 1 + divProtCell + ind] = bdrLoc[numExtraLayers - 2 + divProtCell + ind] + thckInAs/divNCell
	end
	#=
	# Fill boundaries for P-type part of the cell
	for ind = 0 : divPCell - 1
		bdrLoc[4 + divProtCell + divNCell + ind + 1] = bdrLoc[4 + divProtCell + divNCell + ind] + xMinP*(mboxP^ind)
	end
	=#
	for ind = 1 : divPCell

		bdrLoc[numExtraLayers - 1 + divNCell + divProtCell + ind] = bdrLoc[numExtraLayers - 2 + divNCell + divProtCell + ind] + thckInAsSbP/divPCell
	end
	#
	# Fill boundaries for substrate
	for ind = 0 : divSubCell - 1
		bdrLoc[numExtraLayers + divProtCell + divNCell + divPCell + ind ] = bdrLoc[numExtraLayers - 1 + divProtCell + divNCell + divPCell + ind] + xMinSub*(mboxSub^ind)
	end
	#=
	for ind = 1 : divSubCell

		bdrLoc[3 + divNCell + divPCell + ind] = bdrLoc[2 + divNCell + divPCell + ind] + resSubCell
	end
	=#


	## Optical models and absorption functions.
	# Silicon emitter
	if abs(Si_ndoping) == Si_ndoping
		acpDpt = dptDsc(0.045, 0.0)
		dnrDpt = dptDsc(0.0456, abs(Si_ndoping)/1.0e6)
	else
		acpDpt = dptDsc(0.045, abs(Si_ndoping)/1.0e6)
		dnrDpt = dptDsc(0.0456, 0.0)
	end		

	# Calculate silicon model parameters, see siMod and ResponseModels.
	siModParamsE = prmMSi(tEmit, dnrDpt, acpDpt)
	# Construct silicon response model.
	siRspE(enr) = siRsp(enr, siModParamsE)
	siAbsE(enr) = imag(siRspE(enr))

	# Gap response.
	gRsp(enr) = cstRsp(1.0 + im * ^(10.0, -7), enr)
	gAbs(enr) = 0.0


	###abs is the trfFacs to calculate photon numbers, the other imag is to calculate total heat transfer
	#protection layer
	InAsSbPstructure_prot = InP_struct(abs(prot_ndoping),300.0)
	#n-type
	eps_prot_n(enr) = eps_InP_ntype(enr,InAsSbPstructure_prot)
	abs_prot_n_IBimag(enr) = /(eps_InP_imag(enr,InAsSbPstructure_prot), enr)
	eps_prot_n_imag(enr) = imag(eps_InP_ntype(enr,InAsSbPstructure_prot))
	#p-type
	eps_prot_p(enr) = eps_InP_ptype(enr,InAsSbPstructure_prot)
	abs_prot_p_IBimag(enr) = /(eps_InP_imag(enr,InAsSbPstructure_prot), enr)
	eps_prot_p_imag(enr) = imag(eps_InP_ptype(enr,InAsSbPstructure_prot))
	
	
	
	# emitter layer
	InAs_param =InGaAs_struct(abs(ndoping_InAs),300.0)
	eps_emitter_n(enr) = eps_InGaAs_ntype(enr,InAs_param)
	abs_emitter_n_IBimag(enr) = /(eps_InGaAs_imag(enr,InAs_param), enr)
	eps_emitter_n_imag(enr) = imag(eps_InGaAs_ntype(enr,InAs_param))
	#p-type
	eps_emitter_p(enr) = eps_InGaAs_ptype(enr,InAs_param)
	abs_emitter_p_IBimag(enr) = /(eps_InGaAs_imag(enr,InAs_param), enr)
	eps_emitter_p_imag(enr) = imag(eps_InGaAs_ptype(enr,InAs_param))

	

	# base layer
	InAsSbPstructure = InP_struct(abs(pdoping_InAs),300.0)
	#n-type
	eps_base_n(enr) = eps_InP_ntype(enr,InAsSbPstructure)
	abs_base_n_IBimag(enr) = /(eps_InP_imag(enr,InAsSbPstructure), enr)
	eps_base_n_imag(enr) = imag(eps_InP_ntype(enr,InAsSbPstructure))
	#p-type
	eps_base_p(enr) = eps_InP_ptype(enr,InAsSbPstructure)
	abs_base_p_IBimag(enr) = /(eps_InP_imag(enr,InAsSbPstructure), enr)
	eps_base_p_imag(enr) = imag(eps_InP_ptype(enr,InAsSbPstructure))
	

	# InAs p-type subtrate doping and temperature
	pInAs_param = InGaAs_struct(abs(substrate_doping),300.0)
	#n-type
	eps_sub_n(enr) = eps_InGaAs_ntype(enr,pInAs_param)
	abs_sub_n_IBimag(enr) = /(eps_InGaAs_imag(enr,pInAs_param), enr)
	eps_sub_n_imag(enr) = imag(eps_InGaAs_ntype(enr,pInAs_param))
	#p-type
	eps_sub_p(enr) = eps_InGaAs_ptype(enr,pInAs_param)
	abs_sub_p_IBimag(enr) = /(eps_InGaAs_imag(enr,pInAs_param), enr)
	eps_sub_p_imag(enr) = imag(eps_InGaAs_ptype(enr,pInAs_param))
	
	#gold layer epsilon
	eps_gold_imag(enr) = imag(epsgold(enr))


	## Generate lists of optical responses and transfer factors.
	optRsp = []
	trfFacs = []  #absorption

	# Layers prior to N-type part of the cell.
	if firstGap == 0.0
		push!(optRsp, siRspE, gRsp) 
		push!(trfFacs, siAbsE, gAbs)
	else
		push!(optRsp, gRsp, siRspE, gRsp) 
		push!(trfFacs, gAbs, siAbsE, gAbs)
	end
	# Protection layer part of the PV cell
	for ind = 1 : divProtCell

		
		if abs(prot_ndoping) == prot_ndoping 
			#checks doping type
			push!(optRsp,eps_prot_n)
			if Heat_or_Photon == 0
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_prot_n_imag)
			else
				push!(trfFacs,abs_prot_n_IBimag)
			end
		else
			push!(optRsp,eps_prot_p)
			if Heat_or_Photon == 0
				push!(trfFacs,eps_prot_p_imag)
			else
				push!(trfFacs,abs_prot_p_IBimag)
			end
		end
	end

	# N-type part of the PV cell. n-InAs
	for ind = 1 : divNCell

		if abs(ndoping_InAs) == ndoping_InAs
			push!(optRsp, eps_emitter_n)
			if Heat_or_Photon == 0
				push!(trfFacs, eps_emitter_n_imag)
			else
				push!(trfFacs, abs_emitter_n_IBimag)
			end
		else
			push!(optRsp, eps_emitter_p)
			if Heat_or_Photon == 0
				push!(trfFacs, eps_emitter_p_imag)
			else
				push!(trfFacs, abs_emitter_p_IBimag)
			end
		end
	end
	# P-type part of the PV cell. InAsSbP
	for ind = 1 : divPCell
		
		if abs(pdoping_InAs) == pdoping_InAs
			push!(optRsp, eps_base_n)
			if Heat_or_Photon == 0
				push!(trfFacs, eps_base_n_imag)
			else
				push!(trfFacs, abs_base_n_IBimag)
			end
		else
			push!(optRsp, eps_base_p)
			if Heat_or_Photon == 0
				push!(trfFacs, eps_base_p_imag)
			else
				push!(trfFacs, abs_base_p_IBimag)
			end
	
		end
	end
	#P-type substrate of PV cell. InAs
	for ind = 1 : divSubCell

		if abs(substrate_doping) == substrate_doping
			push!(optRsp,eps_sub_n)
			if Heat_or_Photon == 0
				push!(trfFacs,eps_sub_n_imag)
			else
				push!(trfFacs, abs_sub_n_IBimag)
			end
		else
			push!(optRsp,eps_sub_p)
			if Heat_or_Photon == 0
				push!(trfFacs,eps_sub_p_imag)
			else
				push!(trfFacs, abs_sub_p_IBimag)
			end
		end
	end
	# Optically thick backing, gold/p-InAs
	push!(optRsp, epsgold)
	push!(trfFacs,eps_gold_imag)
	## Set which layers transmission should be calculated for.
	if Heat_or_Photon == 0
		NumberlPairs = divProtCell + divNCell + divPCell + divSubCell + 1
	else
		NumberlPairs = divProtCell + divNCell + divPCell + divSubCell
	end
	lPairs = Array{Int64,2}(undef, 2, NumberlPairs)
	#  layer pairs.
	for ind = 1 : NumberlPairs

		lPairs[1, ind] = numExtraLayers - 2
		lPairs[2, ind] =  numExtraLayers - 1 + ind
	end

	# Build layer description for heat transfer code.
	return (lyrDsc(bdrLoc, tmpLst, optRsp, trfFacs),lPairs)
end
precompile(uOttawaSlabs_v3_InGaAs_InP, (Float64, Float64, Int, Int, Int, Int,Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int))

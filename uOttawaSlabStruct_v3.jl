# Construction of PV structure for University of Ottawa collaboration.
# See slides for schematic, slab structure is Si,Gap, Si, Gap, InAsSbP, InAs (N cell),
# InAsSbP (P cell), InAs, gold.
# All temperatures are in Kelvin, and all lengths are in microns.
# Absorption functions are current setup to calculate the number of photons absorbed as interband transitions.
# - v3: counts above interband photon absorption numbers and heat transfer, #Input: Heat_or_Photon , 0 == Heat , 1 == Photon
####


const e = 1.602*^(10,-19)
const eV = 1.60218*^(10,-19)
const hbar = 1.05457*^(10,-34) #m2kg/s
const hbEV = 6.5821*^(10,-16) #eV s
function uOttawaSlabs_v3(tEmit::Float64, tBck::Float64, divProtCell::Int, divNCell::Int, divPCell::Int, divBSFCell::Int, divSubCell::Int, firstGap::Float64, thckRad::Float64, distGap::Float64, thckProt::Float64, thckInAs::Float64, thckInAsSbP::Float64, thckBSF::Float64, thckSub::Float64, Si_ndoping::Float64, prot_ndoping::Float64, ndoping_InAs::Float64, pdoping_InAs::Float64, doping_BSF::Float64, substrate_doping::Float64, quat_x::Float64, prot_quat_x::Float64, quat_BSF_x::Float64, mboxProt::Float64, mboxN::Float64, mboxP::Float64, mboxSub::Float64,xMinProt::Float64, xMinN::Float64, xMinP::Float64, xMinSub::Float64, Heat_or_Photon::Int)
### Settings
	# Total number of layers. + divide radiator to no more than 10 um thick slices
	if firstGap == 0.0
		numExtraLayers = 3
		divRad = 1
	else
		numExtraLayers = 5
		divRad =1#Int(ceil(thckRad/10.0))
	end
	numExtraLayers += divRad -1
	numLayers = divProtCell + divNCell + divPCell + divBSFCell + divSubCell + numExtraLayers
	## Temperature list of the layers.
	tmpLst = fill(tBck, numLayers)
	# Set temperature of the emitter.
	#tmpLst[numExtraLayers-2] = tEmit
	tmpLst[numExtraLayers-1-divRad:numExtraLayers-2] .= tEmit
	## Boundary locations.
	bdrLoc = Vector{Float64}(undef, numLayers - 1) 
	bdrLoc[1] = 0.0
	if firstGap == 0.0
		bdrLoc[2] = distGap
	else
		bdrLoc[2] = firstGap
		for ind = 1 : divRad
			bdrLoc[2 + ind] = bdrLoc[1+ ind] + thckRad/divRad  
		end
		bdrLoc[numExtraLayers-1] = distGap + thckRad + firstGap
	end
	#bdrLoc[5] = bdrLoc[4] + thckProt 
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

	# Fill boundaries of BSF layer 
	for ind = 1 : divBSFCell

		bdrLoc[numExtraLayers - 1 + divPCell + divNCell + divProtCell + ind] = bdrLoc[numExtraLayers - 2 + divNCell + divPCell + divProtCell + ind] + thckBSF/divBSFCell
	end
	#
	# Fill boundaries for substrate
	for ind = 0 : divSubCell - 1
		bdrLoc[numExtraLayers + divBSFCell + divNCell + divPCell + divProtCell + ind] = bdrLoc[numExtraLayers - 1 + divBSFCell + divNCell + divPCell + divProtCell + ind] + xMinSub*(mboxSub^ind)
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
	Dpt = dptDsc(0.0456, 1.0e13) #for no doping

	#calculates silicon model parameters for slab way above
	siModParamsE_top = prmMSi(300.0, Dpt, Dpt)
	# Construct silicon response model.
	siRspE_top(enr) = siRsp(enr, siModParamsE_top)
	siAbsE_top(enr) = imag(siRspE_top(enr))

	# Calculate silicon model parameters, see siMod and ResponseModels.
	siModParamsE = prmMSi(tEmit, dnrDpt, acpDpt)
	# Construct silicon response model.
	siRspE(enr) = siRsp(enr, siModParamsE)
	siAbsE(enr) = imag(siRspE(enr))

	# Gap response.
	gRsp(enr) = cstRsp(1.0 + im * ^(10.0, -7), enr)
	gAbs(enr) = 0.0

	#protection layer
	InAsSbPstructure_prot = InAsSbP_struct(prot_quat_x,0.311*(1 - prot_quat_x),abs(prot_ndoping),300.0)
	#n-type
	eps_prot_n(enr) = eps_InAsSbP_xy_ntype(enr,InAsSbPstructure_prot)
	abs_prot_n_IBimag(enr) = /(eps_InAsSbP_imag_xy_ntype(enr,InAsSbPstructure_prot), enr)
	eps_prot_n_imag(enr) = imag(eps_InAsSbP_xy_ntype(enr,InAsSbPstructure_prot))
	#p-type
	eps_prot_p(enr) = eps_InAsSbP_xy_ptype(enr,InAsSbPstructure_prot)
	abs_prot_p_IBimag(enr) = /(eps_InAsSbP_imag_xy_ptype(enr,InAsSbPstructure_prot), enr)
	eps_prot_p_imag(enr) = imag(eps_InAsSbP_xy_ptype(enr,InAsSbPstructure_prot))
	
	
	
	# emitter layer
	InAs_param =eps_InAs_struct(abs(ndoping_InAs),300.0)
	eps_emitter_n(enr) = eps_InAsntype(enr,InAs_param)
	abs_emitter_n_IBimag(enr) = /(imag(epsIBEV(enr,InAs_param.N0,InAs_param.T,InAs_param.E0_T_InAs_value,InAs_param.eps_inf,InAs_param.P,InAs_param.mstar_ptype_hh,InAs_param.mstar_ptype_lh,InAs_param.F,0.0)), enr)
	eps_emitter_n_imag(enr) = imag(eps_InAsntype(enr,InAs_param))
	#p-type
	eps_emitter_p(enr) = eps_InAsptype(enr,InAs_param)
	abs_emitter_p_IBimag(enr) = /(imag(epsIBEV(enr,InAs_param.N0,InAs_param.T,InAs_param.E0_T_InAs_value,InAs_param.eps_inf,InAs_param.P,InAs_param.mstar_ptype_hh,InAs_param.mstar_ptype_lh,InAs_param.F,1.0)), enr)
	eps_emitter_p_imag(enr) = imag(eps_InAsptype(enr,InAs_param))


	# base layer
	InAsSbPstructure = InAsSbP_struct(quat_x,0.311*(1 - quat_x),abs(pdoping_InAs),300.0)
	#n-type
	eps_base_n(enr) = eps_InAsSbP_xy_ntype(enr,InAsSbPstructure)
	abs_base_n_IBimag(enr) = /(eps_InAsSbP_imag_xy_ntype(enr,InAsSbPstructure), enr)
	eps_base_n_imag(enr) = imag(eps_InAsSbP_xy_ntype(enr,InAsSbPstructure))
	#p-type
	eps_base_p(enr) = eps_InAsSbP_xy_ptype(enr,InAsSbPstructure)
	abs_base_p_IBimag(enr) = /(eps_InAsSbP_imag_xy_ptype(enr,InAsSbPstructure), enr)
	eps_base_p_imag(enr) = imag(eps_InAsSbP_xy_ptype(enr,InAsSbPstructure))
	

	# bsf layer
	InAsSbPstructure_BSF = InAsSbP_struct(quat_BSF_x,0.311*(1 - quat_BSF_x),abs(doping_BSF),300.0)
	#n-type
	eps_bsf_n(enr) = eps_InAsSbP_xy_ntype(enr,InAsSbPstructure_BSF)
	abs_bsf_n_IBimag(enr) = /(eps_InAsSbP_imag_xy_ntype(enr,InAsSbPstructure_BSF), enr)
	eps_bsf_n_imag(enr) = imag(eps_InAsSbP_xy_ntype(enr,InAsSbPstructure_BSF))
	#p-type
	eps_bsf_p(enr) = eps_InAsSbP_xy_ptype(enr,InAsSbPstructure_BSF)
	abs_bsf_p_IBimag(enr) = /(eps_InAsSbP_imag_xy_ptype(enr,InAsSbPstructure_BSF), enr)
	eps_bsf_p_imag(enr) = imag(eps_InAsSbP_xy_ptype(enr,InAsSbPstructure_BSF))

	eps_gold_imag(enr) = imag(epsgold(enr))
	# InAs subtrate doping and temperature
	pInAs_param = eps_InAs_struct(abs(substrate_doping),300.0)
	#n-type
	eps_sub_n(enr) = eps_InAsntype(enr,pInAs_param)
	abs_sub_n_IBimag(enr) = /(imag(epsIBEV(enr,pInAs_param.N0,pInAs_param.T,pInAs_param.E0_T_InAs_value_ptype,pInAs_param.eps_inf,pInAs_param.P,pInAs_param.mstar_ptype_hh,pInAs_param.mstar_ptype_lh,pInAs_param.F,0.0)), enr)
	eps_sub_n_imag(enr) = imag(eps_InAsntype(enr,pInAs_param))
	#p-type
	eps_sub_p(enr) = eps_InAsptype(enr,pInAs_param)
	abs_sub_p_IBimag(enr) = /(imag(epsIBEV(enr,pInAs_param.N0,pInAs_param.T,pInAs_param.E0_T_InAs_value_ptype,pInAs_param.eps_inf,pInAs_param.P,pInAs_param.mstar_ptype_hh,pInAs_param.mstar_ptype_lh,pInAs_param.F,1.0)), enr)
	eps_sub_p_imag(enr) = imag(eps_InAsptype(enr,pInAs_param))


	## Generate lists of optical responses and transfer factors.
	optRsp = []
	trfFacs = []  #absorption

	# Layers prior to the cell.
	if firstGap == 0.0
		push!(optRsp, siRspE, gRsp) 
		push!(trfFacs, siAbsE, gAbs)
	else
		#push!(optRsp, siRspE_top, gRsp) 
		#push!(trfFacs, siAbsE_top, gAbs)
		push!(optRsp, gRsp, gRsp) 
		push!(trfFacs, gAbs, gAbs)
		for ind = 1 : divRad
			push!(optRsp, siRspE) 
			push!(trfFacs, siAbsE)
		end
		push!(optRsp, gRsp) 
		push!(trfFacs, gAbs)
	end
	# FSF layer of the PV cell. InAsSbP
	for ind = 1 : divProtCell

		if abs(prot_ndoping) == prot_ndoping 
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
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_prot_p_imag)
			else
				push!(trfFacs,abs_prot_p_IBimag)
			end
		end
	end

	# Emitter layer of the PV cell. InAs
	for ind = 1 : divNCell

		if abs(ndoping_InAs) == ndoping_InAs
			push!(optRsp, eps_emitter_n)
			if Heat_or_Photon == 0
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_emitter_n_imag)
			else
				push!(trfFacs, abs_emitter_n_IBimag)
			end
		else
			push!(optRsp, eps_emitter_p)
			if Heat_or_Photon == 0
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_emitter_p_imag)
			else
				push!(trfFacs, abs_emitter_p_IBimag)
			end
		end
	end
	# base layer of PV cell. InAsSbP
	for ind = 1 : divPCell
		
		if abs(pdoping_InAs) == pdoping_InAs
			push!(optRsp, eps_base_n)
			if Heat_or_Photon == 0
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_base_n_imag)
			else
				push!(trfFacs, abs_base_n_IBimag)
			end
		else
			push!(optRsp, eps_base_p)
			if Heat_or_Photon == 0
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_base_p_imag)
			else
				push!(trfFacs, abs_base_p_IBimag)
			end
		end
	end
	# bsf layer of PV cell. InAsSbP
	for ind = 1 : divBSFCell
		
		if abs(doping_BSF) == doping_BSF
			push!(optRsp, eps_bsf_n)
			if Heat_or_Photon == 0
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_bsf_n_imag)
			else
				push!(trfFacs, abs_bsf_n_IBimag)
			end
		else
			push!(optRsp, eps_bsf_p)
			if Heat_or_Photon == 0
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_bsf_p_imag)
			else
				push!(trfFacs, abs_bsf_p_IBimag)
			end
		end
	end
	# substrate of PV cell, InAs. Also add the backing layer.
	for ind = 1 : divSubCell + 1

		if ind == divSubCell + 1 && thckSub < 5000.0
			# add gold backing if thin substrate
			push!(optRsp, epsgold)
			push!(trfFacs,eps_gold_imag)
			continue
		end

		if abs(substrate_doping) == substrate_doping
			push!(optRsp,eps_sub_n)
			if Heat_or_Photon == 0
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_sub_n_imag)
			else
				push!(trfFacs, abs_sub_n_IBimag)
			end
		else
			push!(optRsp,eps_sub_p)
			if Heat_or_Photon == 0
				#checks heat transfer or photon count calculation  
				push!(trfFacs,eps_sub_p_imag)
			else
				push!(trfFacs, abs_sub_p_IBimag)
			end
		end
	end

	## Set which layers transmission should be calculated for
	if Heat_or_Photon == 0 || thckSub >= 5000.0
		NumberlPairs = divProtCell + divNCell + divPCell + divBSFCell + divSubCell + 1
	else 
		NumberlPairs = divProtCell + divNCell + divPCell + divBSFCell + divSubCell 
	end
	lPairs = Array{Int64,2}(undef, 2, NumberlPairs*divRad)
	#  layer pairs.
	for ind = 1 : NumberlPairs

		for dRad = 0 : divRad -1 
			lPairs[1, ind + dRad*NumberlPairs] = numExtraLayers - 2 - dRad
			lPairs[2, ind + dRad*NumberlPairs] = numExtraLayers - 1 + ind
		end
	end

	# Build layer description for heat transfer code.
	return (lyrDsc(bdrLoc, tmpLst, optRsp, trfFacs),lPairs)
end
precompile(uOttawaSlabs_v3, (Float64, Float64, Int, Int, Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int))

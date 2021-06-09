# Construction of PV structure for University of Ottawa collaboration.
# See slides for schematic, slab structure is Si,Gap, Si, Gap, InAsSbP, InAs (N cell),
# InAsSbP (P cell), InAs, gold.
# All temperatures are in Kelvin, and all lengths are in microns.
# Absorption functions are current setup to calculate the number of photons absorbed as interband transitions.
# - v3: doping for FSF layer, counts above interband photon absorption numbers 
####


const e = 1.602*^(10,-19)
const eV = 1.60218*^(10,-19)
const hbar = 1.05457*^(10,-34) #m2kg/s
const hbEV = 6.5821*^(10,-16) #eV s
function uOttawaSlabs_v3(tEmit::Float64, tBck::Float64, divProtCell::Int, divNCell::Int, divPCell::Int, divSubCell::Int, firstGap::Float64, thckRad::Float64, distGap::Float64, thckProt::Float64, thckInAs::Float64, thckInAsSbP::Float64, thckSub::Float64, Si_ndoping::Float64, prot_ndoping::Float64, ndoping_InAs::Float64, pdoping_InAs::Float64, quat_x::Float64, prot_quat_x::Float64, mboxProt::Float64, mboxN::Float64, mboxP::Float64, mboxSub::Float64,xMinProt::Float64, xMinN::Float64, xMinP::Float64, xMinSub::Float64)
### Settings
	# Total number of layers.
	numLayers = 5 + divProtCell + divNCell + divPCell + divSubCell
	## Temperature list of the layers.
	tmpLst = fill(tBck, numLayers)
	# Set temperature of the emitter.
	tmpLst[3] = tEmit
	## Boundary locations.
	bdrLoc = Vector{Float64}(undef, numLayers - 1) 
	bdrLoc[1] = 0.0
	bdrLoc[2] = firstGap
	bdrLoc[3] = thckRad +firstGap 
	bdrLoc[4] = distGap + thckRad + firstGap
	#bdrLoc[5] = bdrLoc[4] + thckProt 
	#fill oundaries with protection layer, using the graded layer thcikness scheme
	for ind = 1 : divProtCell
		bdrLoc[4 + ind] = bdrLoc[3  + ind] + thckProt/divProtCell
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

		bdrLoc[4 + divProtCell + ind] = bdrLoc[3 + divProtCell + ind] + thckInAs/divNCell
	end
	#=
	# Fill boundaries for P-type part of the cell
	for ind = 0 : divPCell - 1
		bdrLoc[4 + divProtCell + divNCell + ind + 1] = bdrLoc[4 + divProtCell + divNCell + ind] + xMinP*(mboxP^ind)
	end
	=#
	for ind = 1 : divPCell

		bdrLoc[4 + divNCell + divProtCell + ind] = bdrLoc[3 + divNCell + divProtCell + ind] + thckInAsSbP/divPCell
	end
	#
	# Fill boundaries for substrate
	for ind = 0 : divSubCell - 1
		bdrLoc[4 + divProtCell + divNCell + divPCell + ind + 1] = bdrLoc[4 + divProtCell + divNCell + divPCell + ind] + xMinSub*(mboxSub^ind)
	end
	#=
	for ind = 1 : divSubCell

		bdrLoc[3 + divNCell + divPCell + ind] = bdrLoc[2 + divNCell + divPCell + ind] + resSubCell
	end
	=#


	## Optical models and absorption functions.
	# Silicon emitter
	acpDpt = dptDsc(0.045, 0.0)
	dnrDpt = dptDsc(0.0456, Si_ndoping/1.0e6)  #converting from m-3 to cm-3 for model convention

	#calculates silicon model parameters for slab way above
	siModParamsE_top = prmMSi(300.0, acpDpt, acpDpt)
	# Construct silicon response model.
	siRspE_top(enr) = siRsp(enr, siModParamsE_top)
	siAbsE_top(enr) = imag(siRspE(enr))

	# Calculate silicon model parameters, see siMod and ResponseModels.
	siModParamsE = prmMSi(tEmit, dnrDpt, acpDpt)
	# Construct silicon response model.
	siRspE(enr) = siRsp(enr, siModParamsE)
	siAbsE(enr) = imag(siRspE(enr))

	# Gap response.
	gRsp(enr) = cstRsp(1.0 + im * ^(10.0, -7), enr)
	gAbs(enr) = 0.0

	#protection layer
	# InAsSbP 
	InAsSbPstructure_prot = InAsSbP_struct(prot_quat_x,0.311*(1 - prot_quat_x),prot_ndoping,300.0)
	eps_InAsSbP_prot(enr) = eps_InAsSbP_xy_ntype(enr,InAsSbPstructure_prot)
	abs_InAsSbP_IBimag_prot(enr) = /(eps_InAsSbP_imag_xy_ntype(enr,InAsSbPstructure_prot), enr)
	# InAs
	InAs_param_prot =eps_InAs_struct(prot_ndoping,300.0)
	eps_InAs_prot(enr) = eps_InAsntype(enr,InAs_param)
	abs_InAs_IBimag_prot(enr) = /(imag(epsIBEV(enr,InAs_param_prot.N0,InAs_param_prot.T,InAs_param_prot.E0_T_InAs_value,InAs_param_prot.eps_inf,InAs_param_prot.P,InAs_param_prot.mstar_ptype_hh,InAs_param_prot.mstar_ptype_lh,InAs_param_prot.F,0.0)), enr)

	
	# InAs N-type (PV Cell) doping concentration and temperature
	InAs_param =eps_InAs_struct(ndoping_InAs,300.0)
	eps_nInAs(enr) = eps_InAsntype(enr,InAs_param)
	abs_nInAs_IBimag(enr) = /(imag(epsIBEV(enr,InAs_param.N0,InAs_param.T,InAs_param.E0_T_InAs_value,InAs_param.eps_inf,InAs_param.P,InAs_param.mstar_ptype_hh,InAs_param.mstar_ptype_lh,InAs_param.F,0.0)), enr)
	
	# InAsSbP
	InAsSbPstructure = InAsSbP_struct(quat_x,0.311*(1 - quat_x),pdoping_InAs,300.0)
	eps_InAsSbP(enr) = eps_InAsSbP_xy_ptype(enr,InAsSbPstructure)
	abs_InAsSbP_IBimag(enr) = /(eps_InAsSbP_imag_xy_ptype(enr,InAsSbPstructure), enr)


	eps_gold_imag(enr) = imag(epsgold(enr))
	# InAs p-type subtrate doping and temperature
	pInAs_param = eps_InAs_struct(pdoping_InAs,300.0)
	eps_pInAs(enr) = eps_InAsptype(enr,pInAs_param)
	abs_pInAs_IBimag(enr) = /(imag(epsIBEV(enr,pInAs_param.N0,pInAs_param.T,pInAs_param.E0_T_InAs_value_ptype,pInAs_param.eps_inf,pInAs_param.P,pInAs_param.mstar_ptype_hh,pInAs_param.mstar_ptype_lh,pInAs_param.F,1.0)), enr)
	
	## Generate lists of optical responses and transfer factors.
	optRsp = []
	trfFacs = []  #absorption

	# Layers prior to N-type part of the cell.
	#push!(optRsp, siRspE, gRsp,eps_InAsSbP) 
	push!(optRsp, siRspE_top, gRsp, siRspE, gRsp) 
	#Si, air and protective InAsSbP
	#push!(trfFacs, siAbsE, gAbs,eps_InAsSbP_IBimag)
	push!(trfFacs, siAbsE_top, gRsp, siAbsE, gAbs)
	# Protection layer part of the PV cell
	if prot_quat_x != 1.0
		for ind = 1 : divProtCell

			push!(optRsp,eps_InAsSbP_prot)
			push!(trfFacs,abs_InAsSbP_IBimag_prot)
		end
	else
		for ind = 1 : divProtCell

			push!(optRsp,eps_InAs_prot)
			push!(trfFacs,abs_InAs_IBimag_prot)
		end
	end
	# N-type part of the PV cell. n-InAs
	for ind = 1 : divNCell

		push!(optRsp, eps_nInAs)
		push!(trfFacs, abs_nInAs_IBimag)
	end
	# P-type part of the PV cell. InAsSbP
	for ind = 1 : divPCell

		push!(optRsp, eps_InAsSbP)
		push!(trfFacs, abs_InAsSbP_IBimag)
	end
	#P-type substrate of PV cell. InAs
	for ind = 1 : divSubCell

		push!(optRsp,eps_pInAs)
		push!(trfFacs, abs_pInAs_IBimag)
	end
	# Optically thick backing, gold/p-InAs
	push!(optRsp, epsgold)
	push!(trfFacs,eps_gold_imag)
	## Set which layers transmission should be calculated for.
	lPairs = Array{Int64,2}(undef, 2, divProtCell + divNCell + divPCell + divSubCell)
	#  layer pairs.
	for ind = 1 : divNCell + divPCell + divSubCell + divProtCell

		lPairs[1, ind] = 3
		lPairs[2, ind] = 4 + ind
	end

	# Build layer description for heat transfer code.
	return (lyrDsc(bdrLoc, tmpLst, optRsp, trfFacs),lPairs)
end
precompile(uOttawaSlabs_v3, (Float64, Float64, Int, Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64))

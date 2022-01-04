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
function uOttawaSlabs_v3(layer_temp::Array{Float64,1}, layer_disc::Array{Int64,1} layer_d::Array{Float64,1},layer_dop::Array{Float64,1},layer_Mat::Any,Output::Int64)
### Settings
## 		Input: 	layer_temp: layer temperature (K)
#				layer_disc: number of layer discretizations
#				layer_d:   layer thickness (um)
#				layer_dop: layer doping (negative is p-type)
#				layer_Mat: layer material (string)
#				Output:  1 = heat transfer, 2 = photon #
	# Total number of layers.
	numLayers = sum(layer_disc)
	## Temperature list of the layers.
	tmpLst = layer_temp
	## Boundary locations.
	bdrLoc = Vector{Float64}(undef, numLayers - 1) 
	bdrLoc[1] = 0.0
	for ind = 2 : length(layer_d)
		for ind_layer = 1 : numLayers
			bdrLoc[ind] = bdrLoc[ind-1] + layer_d[ind]/ind_layer
		end
	end



	## Optical models and absorption functions.
	# Silicon emitter
	if abs
		Dpt = dptDsc(0.0456, 1.0e13) #for no doping
	if abs(layer_dop[1]) == layer_dop[1]
		acpDpt = dptDsc(0.045, layer_dop[1]/1.0e6)
		dnrDpt = dptDsc(0.0456, 0.0)  #converting from m-3 to cm-3 for model convention
	else
	if abs
		Dpt = dptDsc(0.0456, 1.0e13) #for no doping

	#calculates silicon model parameters for slab way above
	siModParamsE_top = prmMSi(300.0, Dpt, Dpt)
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
	pInAs_param = eps_InAs_struct(substrate_doping,300.0)
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
precompile(uOttawaSlabs_v3, (Float64, Float64, Int, Int, Int, Int,Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64))

# Construction of PV structure for University of Ottawa collaboration.
# See slides for schematic, slab structure is Si, Gap, InAsSbP, InAs (N cell),
# InAsSbP (P cell), InAs.
# All temperatures are in Kelvin, and all lengths are in microns.
# Absorption functions are current setup to calculate the number of photons
# absorbed as interband transitions.
# Alter code under the ***Optical models and absorption functions*** heading
# for alternative behavior.
####
# Number of division of the N-type InAs layer.
# Resolution of photon generation is thckInAs / divNCell.
# divNCell = 3
# Number of division of the P-type InAsSbP layer.
# divPCell = 3
# Temperature of the emitter.
# tEmit = 800.0
# Background temperature.
# tBck = 300.0

const e = 1.602*^(10,-19)
const eV = 1.60218*^(10,-19)
const hbar = 1.05457*^(10,-34) #m2kg/s
const hbEV = 6.5821*^(10,-16) #eV s
function uOttawaSlabs_v2(tEmit::Float64, tBck::Float64, divNCell::Int, divPCell::Int, divSubCell::Int, distGap::Float64, thckProt::Float64, thckInAs::Float64, thckInAsSbP::Float64, thckSub::Float64, ndoping_InAs::Float64, pdoping_InAs::Float64, quat_x::Float64, prot_quat_x::Float64, mboxN::Float64, mboxP::Float64, mboxSub::Float64, xMinN::Float64, xMinP::Float64, xMinSub::Float64)
	### Settings
	# Separation distances and thickness in uma
	#distGap = 0.1
	# Total InAs slab.
	#thckInAs = 0.5
	# Total InAsSbP slab.
	#thckInAsSbP = 0.5
	# Quaternary thickness. thckQuat = 0.025
	# Protective layer thickness.
	#thckProt = 0.015
	# substrate thickness
	#thckSub = 120
	#divSubCell = 6
	### Derived parameters
	## Layer numbers.
	resNCell = /(thckInAs, divNCell) #divNCell is the number
	resPCell = /(thckInAsSbP, divPCell)
	resSubCell = /(thckSub, divSubCell)
	# Total number of layers.
	numLayers = 4 + divNCell + divPCell + divSubCell
	## Temperature list of the layers.
	tmpLst = fill(tBck, numLayers)
	# Set temperature of the emitter.
	tmpLst[1] = tEmit
	## Boundary locations.
	bdrLoc = Vector{Float64}(undef, numLayers - 1)
	bdrLoc[1] = 0.0
	bdrLoc[2] = distGap
	bdrLoc[3] = bdrLoc[2] + thckProt
	# Fill boundaries for N-type part of the cell.
	for ind = 0 : divNCell - 1
		bdrLoc[3 + ind + 1] = bdrLoc[3 + ind] + xMinN*(mboxN^ind)
	end
	#=
	for ind = 1 : divNCell

		bdrLoc[3 + ind] = bdrLoc[2 + ind] + resNCell
	end
	=#
	# Fill boundaries for P-type part of the cell
	for ind = 0 : divPCell - 1
		bdrLoc[3 + divNCell + ind + 1] = bdrLoc[3 + divNCell + ind] + xMinP*(mboxP^ind)
	end
	#=
	for ind = 1 : divPCell

		bdrLoc[3 + divNCell + ind] = bdrLoc[2 + divNCell + ind] + resPCell
	end
	=#
	# Fill boundaries for substrate
	for ind = 0 : divSubCell - 1
		bdrLoc[3 + divNCell + divPCell + ind + 1] = bdrLoc[3 + divNCell + divPCell + ind] + xMinSub*(mboxSub^ind)
	end
	#=
	for ind = 1 : divSubCell

		bdrLoc[3 + divNCell + divPCell + ind] = bdrLoc[2 + divNCell + divPCell + ind] + resSubCell
	end
	=#
	## Optical models and absorption functions.
	# Silicon emitter.
	# Dopant densities, currently empty?  
	# Francoeur set real part of p-doped Si to be same as real part of top layer. P7
	acpDpt = dptDsc(0.0, 0.0)
	dnrDpt = dptDsc(0.0, 0.0)
	# Calculate silicon model parameters, see siMod and ResponseModels.
	siModParamsE = prmMSi(tmpLst[1], dnrDpt, acpDpt)
	# Construct silicon response model.
	siRspE(enr) = siRsp(enr, siModParamsE)
	siAbsE(enr) = imag(siRspE(enr))
	# Gap response.
	gRsp(enr) = cstRsp(1.0 + im * ^(10.0, -7), enr)
	gAbs(enr) = 0.0
	# InAs N-type (PV Cell) doping concentration and temperature
	#InAs_param =eps_InAs_struct(10.0^23,300.0)
	InAs_param =eps_InAs_struct(ndoping_InAs,300.0)
	eps_nInAs(enr) = eps_InAsntype(enr,InAs_param)
	abs_nInAs_IBimag(enr) = /(imag(epsIBEV(enr,InAs_param.N0,InAs_param.T,InAs_param.E0_T_InAs_value)), enr)
	# InAsSbP
	InAsSbPstructure = InAsSbP_struct(quat_x,0.311*(1 - quat_x))
	eps_InAsSbP(enr) = eps_InAsSbP_xy(enr,InAsSbPstructure)
	eps_InAsSbP_IBimag(enr) = eps_InAsSbP_imag_xy(enr,InAsSbPstructure)
	abs_InAsSbP_IBimag(enr) = /(eps_InAsSbP_imag_xy(enr,InAsSbPstructure), enr)
	if prot_quat_x != 1.0
		# InAsSbP protection layer
		InAsSbPstructure_prot = InAsSbP_struct(prot_quat_x,0.311*(1 - prot_quat_x))
		eps_InAsSbP_prot(enr) = eps_InAsSbP_xy(enr,InAsSbPstructure_prot)
		eps_InAsSbP_IBimag_prot(enr) = eps_InAsSbP_imag_xy(enr,InAsSbPstructure_prot)
		abs_InAsSbP_IBimag_prot(enr) = /(eps_InAsSbP_imag_xy(enr,InAsSbPstructure_prot), enr)
	else
		InAs_param_prot =eps_InAs_struct(prot_ndoping,300.0)
		eps_nInAs_prot(enr) = eps_InAsntype(enr,InAs_param_prot)
		abs_nInAs_IBimag_prot(enr) = /(imag(epsIBEV(enr,InAs_param_prot.N0,InAs_param_prot.T,InAs_param_prot.E0_T_InAs_value)), enr)
	end
	eps_gold_imag(enr) = imag(epsgold(enr))
	# InAs p-type subtrate doping and temperature
	pInAs_param = eps_InAs_struct(pdoping_InAs,300.0)
	eps_pInAs(enr) = eps_InAsptype(enr,pInAs_param)
	abs_pInAs_IBimag(enr) = /(imag(epsIBEV(enr,pInAs_param.N0,pInAs_param.T,pInAs_param.E0_T_InAs_value)), enr)
	
	## Generate lists of optical responses and transfer factors.
	optRsp = []
	trfFacs = []  #absorption
	if prot_quat_x != 1.0
		# Layers prior to N-type part of the cell.
		#push!(optRsp, siRspE, gRsp,eps_InAsSbP) 
		push!(optRsp, siRspE, gRsp,eps_InAsSbP_prot) 
		#Si, air and protective InAsSbP
		#push!(trfFacs, siAbsE, gAbs,eps_InAsSbP_IBimag)
		push!(trfFacs, siAbsE, gAbs,eps_InAsSbP_IBimag_prot)
	else
		# Layers prior to N-type part of the cell.
		#push!(optRsp, siRspE, gRsp,eps_InAsSbP) 
		push!(optRsp, siRspE, gRsp,eps_nInAs_prot) 
		#Si, air and protective InAsSbP
		#push!(trfFacs, siAbsE, gAbs,eps_InAsSbP_IBimag)
		push!(trfFacs, siAbsE, gAbs,abs_nInAs_IBimag_prot)
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
	lPairs = Array{Int64,2}(undef, 2, 1 + divNCell + divPCell + divSubCell)
	#  layer pairs.
	for ind = 1 : divNCell + divPCell + divSubCell + 1

		lPairs[1, ind] = 1
		lPairs[2, ind] = 2 + ind
	end
	#= P-type layer pairs.
	for ind = 1:divPCell

		lPairs[1, ind + divNCell] = 1
		lPairs[2, ind + divNCell] = 3 + divNCell + ind
	end
	# substrate layer pairs
	for ind = 1:divSubCell

		lPairs[1, ind + divPCell + divNCell] = 1
		lPairs[2, ind + divPCell + divNCell] = 3 + divNCell + divPCell + ind
	end
	=#
	# Build layer description for heat transfer code.
	return (lyrDsc(bdrLoc, tmpLst, optRsp, trfFacs),lPairs)
end
precompile(uOttawaSlabs_v2, (Float64, Float64, Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64))

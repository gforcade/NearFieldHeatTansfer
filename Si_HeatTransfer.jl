## heat transfer between doped silicon slabs
## This is to validate the optical model  

#directory to code
savFileDir = pwd()#"/home/gforc034/hs_Sean_02102021_v3"

push!(LOAD_PATH,savFileDir)
using Base.Threads, ProgressMeter, Roots, Base, Printf
using FilmDataStructures, FilmUtilities, ResponseModels
const evJ = 1.6021774232052327e-19


# T [K], d [um]
Si1_T=1000.0
gap_d=0.1
Si2_T=300.0



#changing the parameters from sentaurus units to optical model units
#dop in cm-3
Si1_dop = 1.0e20
Si2_dop = 1.0e20



####running simulations
#file name
fName = "SiRad"#"Rad"*string(trunc(Int,round(Rad_T)))*"fstGap"*string(trunc(Int,round(firstgap_d*1000)))*"Si"*string(trunc(Int,round(Rad_d*1000)))*"Dop"*string(round(Rad_dop/1000000,sigdigits=1))*"Gap"*string(trunc(Int,round(gap_d*1000)))*"fsf"*string(trunc(Int,round(fsf_d*1000)))*"As"*string(round(fsf_As,sigdigits=1))*"Dop"*string(round(fsf_dop/1000000,sigdigits=1))*"InAs"*string(trunc(Int,round(emitter_d*1000)))*"Dop"*string(round(emitter_dop/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(base_d*1000)))*"As"*string(round(base_As,sigdigits=1))*"Dop"*string(round(base_dop/1000000,sigdigits=1))*"Sub"*string(trunc(Int,round(substrate_d*1000)))*"Dop"*string(round(substrate_dop/1000000,sigdigits=1))*".txt"
#only simulate if file not present
#if !isfile(savFileDir*"/Photon Number/"*fName)
####simulations for energy transfer

#temperature of layers
numLayers = 3
tmpLst = fill(Si2_T, numLayers)
tmpLst[1] = Si1_T
#border Location
bdrLoc = Vector{Float64}(undef, numLayers - 1) 
bdrLoc[1] = 0.0
bdrLoc[2] = gap_d
#silicon dielectric model parameters
acpDpt_top = dptDsc(0.045, 0.0)
dnrDpt_top = dptDsc(0.0456, abs(Si1_dop))
#for bottom
acpDpt = dptDsc(0.045, 0.0)
dnrDpt = dptDsc(0.0456, abs(Si2_dop))
#calculates silicon model parameters for slab way above
siModParamsE_top = prmMSi(Si1_T, dnrDpt_top, acpDpt_top)
# Construct silicon response model.
siRspE_top(enr) = siRsp(enr, siModParamsE_top)
siAbsE_top(enr) = imag(siRspE_top(enr))
#Si for bottom layer
siModParamsE = prmMSi(Si2_T, dnrDpt, acpDpt)
# Construct silicon response model.
siRspE(enr) = siRsp(enr, siModParamsE)
siAbsE(enr) = imag(siRspE(enr))
# Gap response.
gRsp(enr) = cstRsp(1.0 + im * ^(10.0, -7), enr)
gAbs(enr) = 0.0

## Generate lists of optical responses and transfer factors.
optRsp = []
trfFacs = []  #absorption
push!(optRsp, siRspE_top, gRsp, siRspE) 
#air,Si, air, Si, air 
push!(trfFacs, siAbsE_top, gAbs, siAbsE)

## Set which layers transmission should be calculated for.
lPairs = Array{Int64,2}(undef, 2, 1) #+1 for gold back-reflector
lPairs[1, 1] = 1
lPairs[2, 1] = 3
lVar = lyrDsc(bdrLoc, tmpLst, optRsp, trfFacs)

enrRng = (6.6e-4, 1.319) #energy interval
# Builds slab structure, generating lVar and lPairs.
# Arguments are: temperature of the emitter, background temperature, divisions of NCell



# Storage for heat transfer computations.
htPairs = Array{Float64,1}(undef, size(lPairs)[2])
# Compute heat transfer using heat transfer function.
stats2 = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
# Open output file.
#=
fileStream = open(savFileDir*"/Total Heat Transfer/"*fName,"w")
# Double check thread initialization.
write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
## External program settings, internal setting contained in uOttawaSlabStruct.jl
## Begin program execution.
write(fileStream, "Execution:\n")
write(fileStream, "Layer structure constructed in " * string(stats1.time) * " s.\n")
write(fileStream, "Generated excitation profile computed " * string(stats2.time) * " s.\n")
write(fileStream,string(xMinN)*" "*string(xMinP)*" 0.0 "*string(xMinProt)*" xMin for mbox meshing. Emitter, Base, Substrate, Prot respectively. \n")
write(fileStream,string(divNcell)*" "*string(divPcell)*" "*string(divSubCell)*" "*string(divProtCell)*" of slices for emitter, base, and Substrate, Prot respectively.\n\n")
# Write results to file.
write(fileStream, "Layer boundary edge (microns)	"*" Excitation Density Rate (W cm-2 microns-1)\n\n")
# Location of cell boundary within the slab structure.
# Output calculation results.
for ind = 1 : length(htPairs) - 1
    write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(htPairs[ind]/((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind]-1])),sigdigits=4)) * " \n")
end
#gold backing absorption
write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,length(htPairs)-1]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(htPairs[length(htPairs)],sigdigits=4)) * " \n")
# Flush file stream.
close(fileStream)
=#
println("Heattransfer execution time = ",stats2.time)
println("Total heat transfer (W/m2) = ", @sprintf("%.2E",sum(htPairs)*1e4))
    


   
#end









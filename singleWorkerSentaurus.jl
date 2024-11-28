# - v3: calculates the total depth resolved heat transfer + photon number to cell
# new structure

#directory to code
savFileDir = "/home/gforc034/hs_Sean_02102021_v3"

push!(LOAD_PATH,savFileDir)
using Base.Threads, ProgressMeter, Roots, Base,Printf
using FilmDataStructures, FilmUtilities, ResponseModels, OttawaSlabCell_v3
const evJ = 1.6021774232052327e-19


# T [K], d [um], dop [m-3]
Rad_T=700.0
firstgap_d=0.1
Rad_d=30.0
gap_d=0.1
fsf_d=0.2
emitter_d=3.98
base_d=0.032
bsf_d=5.0
substrate_d=500.0
Rad_dop=-2.0*(10.0^19)
fsf_dop=-1.0*(10.0^19)
emitter_dop=6.0*(10.0^14)
base_dop=1.0*(10.0^16)  
bsf_dop=1e20
substrate_dop=2.0*(10.0^18)
fsf_As=1.0
base_As=1.0  #smallest is best because of higher bandgap
bsf_As=1.0


#changing the parameters from sentaurus units to optical model units
#dop from cm-3 to m-3
Rad_dop = Rad_dop*1e6
fsf_dop = fsf_dop*1e6
emitter_dop = emitter_dop*1e6
base_dop = base_dop*1e6
bsf_dop = bsf_dop*1e6
substrate_dop = substrate_dop*1e6




    
#gets the mbox value for a given N=# fo  slices and r = xthck/xmin
function mboxish(x,N,r)
    
    tot = 0.0
    for ind = 1 : N - 1
        tot += x^ind
    end    
    return tot - (r - 1)
end




####running simulations
#file name
#fName = "Rad"*string(trunc(Int,round(Rad_T)))*"fstGap"*string(trunc(Int,round(firstgap_d*1000)))*"Si"*string(trunc(Int,round(Rad_d*1000)))*"Dop"*string(round(Rad_dop/1000000,sigdigits=1))*"Gap"*string(trunc(Int,round(gap_d*1000)))*"fsf"*string(trunc(Int,round(fsf_d*1000)))*"As"*string(round(fsf_As,sigdigits=1))*"Dop"*string(round(fsf_dop/1000000,sigdigits=1))*"InAs"*string(trunc(Int,round(emitter_d*1000)))*"Dop"*string(round(emitter_dop/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(base_d*1000)))*"As"*string(round(base_As,sigdigits=1))*"Dop"*string(round(base_dop/1000000,sigdigits=1))*"Sub"*string(trunc(Int,round(substrate_d*1000)))*"Dop"*string(round(substrate_dop/1000000,sigdigits=1))*".txt"
#only simulate if file not present

####simulations for energy transfer
xMinProt = 0.0
xMinN = 0.0
xMinP = 0.0
divProtCell = 1
divNcell = 1
divPcell = 1 
if bsf_d != 0.0
    divBSFCell = 1
else
    divBSFCell = 0
end
if substrate_d == 0.0 || substrate_d > 5000.0
    divSubCell = 0 
    xMinSub = 0.0
    mboxSub = 0.0
else
    divSubCell = 1
    xMinSub = substrate_d
    mboxSub = 1.0
end
mboxProt = 0.0
mboxN = 0.0
mboxP = 0.0
enrRng = (0.01, 1.0) #energy interval
# Builds slab structure, generating lVar and lPairs.
# Arguments are: temperature of the emitter, background temperature, divisions of NCell
stats1 = @timed (lVar, lPairs) = uOttawaSlabs_v3(Rad_T, 300.0, divProtCell, divNcell, divPcell, divBSFCell, divSubCell, firstgap_d, Rad_d, gap_d, fsf_d, emitter_d, base_d, bsf_d, substrate_d, Rad_dop, fsf_dop, emitter_dop, base_dop,bsf_dop,substrate_dop, base_As, fsf_As, bsf_As, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub,0)
# Storage for photon number computations.
htPairs = Array{Float64,1}(undef, size(lPairs)[2])
# Compute number of generated photons using heat transfer function.
# Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
stats2 = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
# Open output file.
fileStream = open(savFileDir*"/Total Heat Transfer.txt","w")
# Double check thread initialization.
write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
## External program settings, internal setting contained in uOttawaSlabStruct.jl
## Begin program execution.
write(fileStream, "Execution:\n")
write(fileStream, "Layer structure constructed in " * string(stats1.time) * " s.\n")
write(fileStream, "Generated excitation profile computed " * string(stats2.time) * " s.\n")
write(fileStream,string(xMinN)*" "*string(xMinP)*" "*string(xMinSub)*" "*string(xMinProt)*" xMin for mbox meshing. InAs, Q, Substrate, Prot respectively. \n")
write(fileStream,string(divNcell)*" "*string(divPcell)*" "*string(divSubCell)*" "*string(divProtCell)*" "*string(divBSFCell)*" of slices for InAs, Q, and Substrate, Prot, BSF respectively.\n\n")
# Write results to file.
write(fileStream, "Layer boundary edge (microns)	"*" Excitation Density Rate (W cm-2 )\n\n")
# Location of cell boundary within the slab structure.
numLayers = divNcell + divPcell + divProtCell + divBSFCell + divSubCell
# Output calculation results.
for ind = 1 : numLayers
    write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(sum(htPairs[ind:numLayers+1:length(htPairs)]),sigdigits=4)) * " \n")
end
#gold backing absorption
write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,length(htPairs)-1]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(sum(htPairs[numLayers+1:numLayers+1:length(htPairs)]),sigdigits=4)) * " \n")
# Flush file stream.
close(fileStream)




###### simulation for photon number
##simulation run with no subdivisions for accurate absolute values
mboxSub = 1.0
#energy interval
enrRng = (0.35, 1.0)
# Builds slab structure, generating lVar and lPairs.
# Arguments are: temperature of the emitter, background temperature, divisions of NCell
stats = @timed (lVar, lPairs) = uOttawaSlabs_v3(Rad_T, 300.0, divProtCell, divNcell, divPcell, divBSFCell, divSubCell, firstgap_d, Rad_d, gap_d, fsf_d, emitter_d, base_d, bsf_d, substrate_d, Rad_dop, fsf_dop, emitter_dop, base_dop,bsf_dop,substrate_dop, base_As, fsf_As, bsf_As, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub,1)
# Storage for photon number computations.
htPairs = Array{Float64,1}(undef, size(lPairs)[2])
# Compute number of generated photons using heat transfer function.
# Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
stats = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
#copy ouput into new array to use later
htPairsNorm = copy(htPairs)
numLayersNorm = divNcell + divPcell + divProtCell + divBSFCell + divSubCell


divProtCell = ceil(Int,sqrt(fsf_d*100))
divNcell = ceil(Int,sqrt(emitter_d*100))
divPcell = ceil(Int,sqrt(base_d*100))
divBSFCell = ceil(Int,sqrt(bsf_d*100))
if bsf_dop == abs(bsf_dop) || divBSFCell == 1
    if divBSFCell != 0
        divBSFCell = 2
    end
end
mboxSuber(x) = mboxish(x,divSubCell,substrate_d/xMinSub)
if substrate_d == 0.0
    divSubCell = 0 
    xMinSub = 0.0
    mboxSub = 1.0
elseif substrate_d < 50.0
    divSubCell = ceil(Int,sqrt(substrate_d*100))
    xMinSub = 0.01
    mboxSub = find_zero(mboxSuber,2.0,Order5())
elseif bsf_dop > 1e19*1e6 
    divSubCell = 3
    xMinSub = 10.0
    mboxSub = find_zero(mboxSuber,10.0,Order5())
else
    divSubCell = 30
    xMinSub = 0.1
    mboxSub = find_zero(mboxSuber,2.0,Order5())
end
## External program settings, internal setting contained in uOttawaSlabStruct.jl
enrRng = (0.35, 1.0)
# Builds slab structure, generating lVar and lPairs.
# Arguments are: temperature of the emitter, background temperature, divisions of NCell
stats1 = @timed (lVar, lPairs) = uOttawaSlabs_v3(Rad_T, 300.0, divProtCell, divNcell, divPcell, divBSFCell, divSubCell, firstgap_d, Rad_d, gap_d, fsf_d, emitter_d, base_d, bsf_d, substrate_d, Rad_dop, fsf_dop, emitter_dop, base_dop,bsf_dop,substrate_dop, base_As, fsf_As, bsf_As, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub,1)
# Storage for photon number computations.
htPairs = Array{Float64,1}(undef, size(lPairs)[2])
# Compute number of generated photons using heat transfer function.
# Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
stats2 = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
# Open output file.
fileStream = open(savFileDir*"/Photon Number.txt","w")
#end
# Double check thread initialization.
write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
## Begin program execution.
write(fileStream, "Execution:\n")
write(fileStream, "Layer structure constructed in " * string(stats1.time) * " s.\n")
write(fileStream, "Generated excitation profile computed " * string(stats2.time) * " s.\n")
write(fileStream,string(xMinN)*" "*string(xMinP)*" "*string(xMinSub)*" "*string(xMinProt)*" xMin for mbox meshing. InAs, Q, Substrate, Prot respectively. \n")
write(fileStream,string(divNcell)*" "*string(divPcell)*" "*string(divSubCell)*" "*string(divProtCell)*" "*string(divBSFCell)*" of slices for InAs, Q, and Substrate, Prot, BSF respectively.\n\n")
# Write results to file.
write(fileStream, "Depth (microns)	"*" Excitation Density Rate (cm-2 microns-1 s-1)\n\n")
idxNorm = 1
extremes = [1,divProtCell]
numLayers = divNcell + divPcell + divProtCell + divBSFCell + divSubCell
# Output calculation results.
for ind = 1 : numLayers
        #calculate total absorption within material layer i
        normhtPairs = 0.0
        for x = extremes[1]:extremes[2] 
            normhtPairs += sum(htPairs[x:numLayers:length(htPairs)])
        end
        write(fileStream, string(round(sum(htPairs[ind:numLayers:length(htPairs)])/(normhtPairs*(lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind]-1])*evJ)*sum(htPairsNorm[idxNorm:numLayersNorm:length(htPairsNorm)]),sigdigits=4)) * " \n")
    #to go through all the layers
    if ind == divProtCell || ind == divProtCell + divNcell || ind == divProtCell + divNcell + divPcell
        global idxNorm += 1 
    end
    if ind == divProtCell
        global extremes = [ind+1,ind + divNcell]
    elseif ind == divProtCell + divNcell
        global extremes = [ind+1, ind + divPcell]
    elseif ind == divProtCell + divNcell + divPcell && divBSFCell != 0
        global extremes = [ind+1,ind+divBSFCell]
    elseif ind == divProtCell + divNcell + divPcell + divBSFCell
        global extremes = [ind+1,ind+divSubCell]
    end
end
# Flush file stream.
close(fileStream)




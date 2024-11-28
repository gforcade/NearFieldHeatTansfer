# - v3: calculates the total depth resolved heat transfer + photon number to cell
# new structure

#directory to code
savFileDir = "/home/gforc034/hs_Sean_02102021_v3"

push!(LOAD_PATH,savFileDir)
using Base.Threads, ProgressMeter, Roots, Base, Printf
using FilmDataStructures, FilmUtilities, ResponseModels, OttawaSlabCell_v3_InGaAs_InP
const evJ = 1.6021774232052327e-19


# T [K], d [um], dop [m-3]
Rad_T=930.0
firstgap_d=0.01
Rad_d=60.0
gap_d=0.1
fsf_d=0.2
emitter_d=1.0
base_d=0.1
substrate_d=0.1
Rad_dop=-3.0*(10.0^19)
fsf_dop=-1.0*(10.0^18)
emitter_dop=1.0*(10.0^17)
base_dop=1.0*(10.0^18)  
substrate_dop=1.0*(10.0^18)
fsf_As=0.47
base_As=0.47  


#changing the parameters from sentaurus units to optical model units
#dop from cm-3 to m-3
Rad_dop = Rad_dop*1e6
fsf_dop = fsf_dop*1e6
emitter_dop = emitter_dop*1e6
base_dop = base_dop*1e6
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
fName = "Rad"*string(trunc(Int,round(Rad_T)))*"fstGap"*string(trunc(Int,round(firstgap_d*1000)))*"Si"*string(trunc(Int,round(Rad_d*1000)))*"Dop"*string(round(Rad_dop/1000000,sigdigits=1))*"Gap"*string(trunc(Int,round(gap_d*1000)))*"fsf"*string(trunc(Int,round(fsf_d*1000)))*"As"*string(round(fsf_As,sigdigits=1))*"Dop"*string(round(fsf_dop/1000000,sigdigits=1))*"InAs"*string(trunc(Int,round(emitter_d*1000)))*"Dop"*string(round(emitter_dop/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(base_d*1000)))*"As"*string(round(base_As,sigdigits=1))*"Dop"*string(round(base_dop/1000000,sigdigits=1))*"Sub"*string(trunc(Int,round(substrate_d*1000)))*"Dop"*string(round(substrate_dop/1000000,sigdigits=1))*".txt"
#only simulate if file not present
if !isfile(savFileDir*"/Photon Number_InGaAs/"*fName)
    ####simulations for energy transfer
    xMinProt = 0.0
    xMinN = 0.0
    xMinP = 0.0
    xMinSub = substrate_d
    divProtCell = 1
    divNcell = 1
    divPcell = 1 
    if substrate_d == 0.0
        divSubCell = 0 
    else 
        divSubCell = 1 
    end
    mboxProt = 0.0
    mboxN = 0.0
    mboxP = 0.0
    mboxSub = 1.0
    enrRng = (0.01, 1.5) #energy interval
    
    # Builds slab structure, generating lVar and lPairs.
    # Arguments are: temperature of the emitter, background temperature, divisions of NCell
    stats1 = @timed (lVar, lPairs) = uOttawaSlabs_v3_InGaAs_InP(Rad_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, firstgap_d, Rad_d, gap_d, fsf_d, emitter_d, base_d, substrate_d, Rad_dop, fsf_dop, emitter_dop, base_dop,substrate_dop, base_As, fsf_As, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub,0)
    # Storage for heat transfer computations.
    htPairs = Array{Float64,1}(undef, size(lPairs)[2])
    # Compute heat transfer using heat transfer function.
    stats2 = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
    # Open output file.
    fileStream = open(savFileDir*"/Total Heat Transfer_InGaAs/"*fName,"w")
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
    numLayers = divNcell + divPcell + divProtCell + divSubCell
    # Output calculation results.
    for ind = 1 : numLayers

        write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(sum(htPairs[ind:numLayers+1:length(htPairs)])/((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind]-1])),sigdigits=4)) * " \n")
    end
    #gold backing absorption
    write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,length(htPairs)-1]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(sum(htPairs[numLayers+1:numLayers+1:length(htPairs)]),sigdigits=4)) * " \n")
    # Flush file stream.
    close(fileStream)
    



    ###### simulation for photon number
    ## simulation run with no subdivisions for accurate absolute values
    #energy interval
    enrRng = (0.685, 1.5)
    # Builds slab structure, generating lVar and lPairs.
    # Arguments are: temperature of the emitter, background temperature, divisions of NCell
    stats = @timed (lVar, lPairs) = uOttawaSlabs_v3_InGaAs_InP(Rad_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, firstgap_d, Rad_d, gap_d, fsf_d, emitter_d, base_d, substrate_d, Rad_dop, fsf_dop, emitter_dop, base_dop, substrate_dop, base_As, fsf_As, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub,1)
    # Storage for photon number computations.
    htPairs = Array{Float64,1}(undef, size(lPairs)[2])
    # Compute number of generated photons using heat transfer function.
    # Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
    stats = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
    #copy ouput into new array to use later
    htPairsNorm = copy(htPairs)
    numLayersNorm = divNcell + divPcell + divProtCell + divSubCell


    #total current
    println( "Simulation runtime: " * string(stats.time) * " s.")
    println("Jph (uA) = ", @sprintf("%.2E",sum(htPairsNorm)*1e6*0.0075^2*pi))   #0.0075 is the nearfield device radius 
    
    
    ## simulation run with subdivisions for depth resolved currrent generation
    divProtCell = ceil(Int,sqrt(fsf_d*100))
    divNcell = ceil(Int,sqrt(emitter_d*100))
    divPcell = ceil(Int,sqrt(base_d*100))
    if substrate_d == 0.0
        divSubCell = 0 
        mboxSub = 1.0 
    elseif substrate_d < 5.0
        divSubCell = ceil(Int,sqrt(substrate_d*100))
        mboxSub = 1.0
        xMinSub = substrate_d/divSubCell
    else 
        #run with xbox method
        if substrate_d > 10.0
            xMinSub = 0.1
            divSubCell = 30
        elseif substrate_d > 5.0
            divSubCell = 30
            xMinSub = 0.01
        end
        mboxSuber(x) = mboxish(x,divSubCell,substrate_d/xMinSub)
        mboxSub = find_zero(mboxSuber,2.0,Order5())
    end
    ## External program settings, internal setting contained in uOttawaSlabStruct.jl
    # Builds slab structure, generating lVar and lPairs.
    # Arguments are: temperature of the emitter, background temperature, divisions of NCell
    stats1 = @timed (lVar, lPairs) = uOttawaSlabs_v3_InGaAs_InP(Rad_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, firstgap_d, Rad_d, gap_d, fsf_d, emitter_d, base_d, substrate_d, Rad_dop, fsf_dop, emitter_dop, base_dop, substrate_dop, base_As, fsf_As, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub,1)
    # Storage for photon number computations.
    htPairs = Array{Float64,1}(undef, size(lPairs)[2])
    # Compute number of generated photons using heat transfer function.
    # Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
    stats2 = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
    # Open output file.
    fileStream = open(savFileDir*"/Photon Number_InGaAs/"*fName,"w")
    # Double check thread initialization.
    write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
    ## Begin program execution.
    write(fileStream, "Execution:\n")
    write(fileStream, "Layer structure constructed in " * string(stats1.time) * " s.\n")
    write(fileStream, "Generated excitation profile computed " * string(stats2.time) * " s.\n")
    write(fileStream,string(xMinN)*" "*string(xMinP)*" "*string(xMinSub)*" "*string(xMinProt)*" xMin for mbox meshing. InAs, Q, Substrate, Prot respectively. \n")
    write(fileStream,string(divNcell)*" "*string(divPcell)*" "*string(divSubCell)*" "*string(divProtCell)*" of slices for InAs, Q, and Substrate, Prot respectively.\n\n")
    # Write results to file.
    write(fileStream, "Depth (microns)	"*" Excitation Density Rate (cm-2 microns-1 s-1)\n\n")
    idxNorm = 1
    extremes = [1,divProtCell]
    numLayers = divNcell + divPcell + divProtCell + divSubCell
    # Output calculation results.
    for ind = 1 : numLayers
        #calculate absorption within material layer i
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
        elseif ind == divProtCell + divNcell + divPcell
            global extremes = [ind+1,ind+divSubCell]
        end
    end
    # Flush file stream.
    close(fileStream)
    
end









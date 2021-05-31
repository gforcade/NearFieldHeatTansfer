push!(LOAD_PATH,pwd())

using Distributed
using Base.Threads, ProgressMeter, Plots, LaTeXStrings
using FilmDataStructures, FilmUtilities, ResponseModels, OttawaSlabCell_v2
const evJ = 1.6021774232052327e-19



# (CHANGE) sweeping over 2 parameters. Change parameters in heatLayersFunc as needed
sweep1 =  [0.1] #LinRange(0.5,3.0,6)
sweep2 = [1.0,3.0,100.0,300.0,1000.0]*(10.0^22)


function heatLayersFunc(d_InAs_total::Float64,ndoping_InAs::Float64)
    #function heatLayersFunc(d_InAs_total::Float64)
    
    d_gap = 0.1
    d_prot=0.015
    #d_InAs_total = 0.5
    d_InAsSbP_base = 0.5
    d_InAs_sub = 120.0
    #ndoping_InAs = 1.0*(10.0^23)
    quat_x = 0.6
    prot_quat_x = 0.6
    divNcell = ceil(Int,d_InAs_total/0.1)
    divPcell = ceil(Int,d_InAsSbP_base/0.1)
    divSubCell = 6
    pdoping_InAs = 3.0*(10.0^24)
    
    
    # Open output file.
    #fileStream = open("photonProfile"*string(divNcell)*".txt","w")
    fileStream = open("Results/SiGap"*string(trunc(Int,round(d_gap*1000)))*"fsf"*string(trunc(Int,round(d_prot*1000)))*"As"*string(prot_quat_x)*"InAs"*string(trunc(Int,round(d_InAs_total*1000)))*"Dop"*string(round(ndoping_InAs/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(d_InAsSbP_base*1000)))*"As"*string(quat_x)*"Sub"*string(trunc(Int,round(d_InAs_sub*1000)))*"Dop"*string(round(pdoping_InAs/1000000,sigdigits=1))*".txt","w")
    # Double check thread initialization.
    write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
    # Conversion of photon energy from electron volts to Joules.
    #const evJ = 1.6021774232052327e-19
    ## External program settings, internal setting contained in uOttawaSlabStruct.jl
    enrRng = (0.06, 1.0)
    ## Begin program execution.
    write(fileStream, "Execution:\n")
    # Builds slab structure, generating lVar and lPairs.
    # Arguments are: temperature of the emitter, background temperature, divisions of NCell
    # (Resolution of photon generation is thckInAs / divNCell), division of PCell
    #stats = @timed (lVar, lPairs) = uOttawaSlabs(800.0, 300.0, divNcell, 30)
    stats = @timed (lVar, lPairs) = uOttawaSlabs_v2(800.0, 300.0, divNcell, divPcell, divSubCell, d_gap, d_prot, d_InAs_total, d_InAsSbP_base, d_InAs_sub, ndoping_InAs, pdoping_InAs, quat_x, prot_quat_x)
    write(fileStream, "Layer structure constructed in " * string(stats.time) * " s.\n")
    # Storage for photon number computations.
    htPairs = Array{Float64,1}(undef, size(lPairs)[2])
    # Compute number of generated photons using heat transfer function.
    # Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
    stats = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
    write(fileStream, "Generated excitation profile computed " * string(stats.time) * " s.\n\n")
    write(fileStream,string(divNcell)*" "*string(divPcell)*" "*string(divSubCell)*" of slices for InAs, Q, and Substrate respectively.\n\n")
    # Write results to file.
    write(fileStream, "Depth (microns)	"*" Excitation Density Rate (cm-2 microns-1 s-1)\n\n")
    # Location of cell boundary within the slab structure.
    cellSurfLoc = lVar.bdrLoc[lPairs[2,1]]
    # Output calculation results.
    for ind = 1 : length(htPairs)
        #htPairs[ind] = htPairs[ind]*(lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind-1]])
        #write(fileStream, string(round(lVar.bdrLoc[lPairs[2,ind]] - cellSurfLoc, sigdigits = 4)) * "," * string(round(/(htPairs[ind], evJ), sigdigits = 4)) * "\n")
        write(fileStream, string(round(htPairs[ind]/((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind]-1])*evJ),sigdigits=4)) * " \n")
    end
    # Flush file stream.
    close(fileStream)
end




# DO NOT NEED TO CHANGE

#function to set up the 2 parameter problem as a matrix
function twoParams(i::Int)
    twosweep = ceil(Int,i/length(sweep1))
    heatLayersFunc(sweep1[i - length(sweep1)*(twosweep-1)],sweep2[twosweep])
end


@showprogress map(i -> twoParams(i),1:length(sweep1)*length(sweep2))

# - v3: calculates the total depth resolved heat transfer + photon number to cell
# new structure

@everywhere push!(LOAD_PATH,pwd())
@everywhere begin
        
    using Distributed
    using Base.Threads, ProgressMeter, Roots
    using FilmDataStructures, FilmUtilities, ResponseModels, OttawaSlabCell_v3_1, OttawaSlabCell_v2, OttawaSlabCell_v3, OttawaSlabCell_v2_1
    const evJ = 1.6021774232052327e-19
end




# (CHANGE) sweeping over 2 parameters. Change parameters in heatLayersFunc as needed
@everywhere sweep1=LinRange(5.0,30.0,6) 
@everywhere sweep2=[3.0, 6.0, 10.0, 30.0, 60.0, 100.0, 300.0, 600.0, 1000.0]*(10.0^24) 


@everywhere function heatLayersFunc(var1::Float64,var2::Float64)
    
    #gets the mbox value for a given N=# fo  slices and r = xthck/xmin
    function mboxish(x,N,r)
        
        tot = 0.0
        for ind = 1 : N - 1
            tot += x^ind
        end
        return tot - (r - 1)
    end


    # T [K], d [um], dop [m-3]
    Rad_T=700.0
    firstgap_d=0.5
    Rad_d=var1 #10.0
    gap_d=0.1
    fsf_d=0.1
    emitter_d=5.0
    base_d=2.0
    substrate_d=120.0
    Rad_dop=var2 #1.0*(10.0^25)
    fsf_dop=3.0*(10.0^24)
    emitter_dop=3.0*(10.0^21)
    base_dop=3.0*(10.0^24)  #highest is best (proven)
    substrate_dop=3.0*(10.0^24)
    fsf_As=1.0
    base_As=0.4  #smallest is best because of higher bandgap


    #changing the parameters from sentaurus units to optical model units
    #dop from cm-3 to m-3
    Rad_dop = Rad_dop*1e6
    fsf_dop = fsf_dop*1e6
    emitter_dop = emitter_dop*1e6
    base_dop = base_dop*1e6
    substrate_dop = substrate_dop*1e6

    ####running simulations
    #file name
    fName = "Rad"*string(trunc(Int,round(Rad_T)))*"fstGap"*string(trunc(Int,round(firstgap_d*1000)))*"Si"*string(trunc(Int,round(Rad_d*1000)))*"Dop"*string(round(Rad_dop/1000000,sigdigits=1))*"Gap"*string(trunc(Int,round(gap_d*1000)))*"fsf"*string(trunc(Int,round(fsf_d*1000)))*"As"*string(round(fsf_As,sigdigits=1))*"Dop"*string(round(fsf_dop/1000000,sigdigits=1))*"InAs"*string(trunc(Int,round(emitter_d*1000)))*"Dop"*string(round(emitter_dop/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(base_d*1000)))*"As"*string(round(base_As,sigdigits=1))*"Dop"*string(round(base_dop/1000000,sigdigits=1))*"Sub"*string(trunc(Int,round(substrate_d*1000)))*"Dop"*string(round(substrate_dop/1000000,sigdigits=1))*".txt"
    #only simulate if file not present
    if !isfile("Total Heat Transfer/"*fName)
        ####simulations for energy transfer
        xMinProt = 0.0
        xMinN = 0.0
        xMinP = 0.0
        xMinSub = 0.0
        divProtCell = 1
        divNcell = 1
        divPcell = 1 
        divSubCell = 1 
        mboxProt = 0.0
        mboxN = 0.0
        mboxP = 0.0
        mboxSub = 0.0
        # Open output file.
        fileStream = open("Total Heat Transfer/"*fName,"w")
        # Double check thread initialization.
        write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
        ## External program settings, internal setting contained in uOttawaSlabStruct.jl
        enrRng = (0.01, 1.0)
        ## Begin program execution.
        write(fileStream, "Execution:\n")
        # Builds slab structure, generating lVar and lPairs.
        # Arguments are: temperature of the emitter, background temperature, divisions of NCell
        stats = @timed (lVar, lPairs) = uOttawaSlabs_v3_1(Rad_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, firstgap_d, Rad_d, gap_d, fsf_d, emitter_d, base_d, substrate_d, Rad_dop, fsf_dop, emitter_dop, base_dop, base_As, fsf_As, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub)
        write(fileStream, "Layer structure constructed in " * string(stats.time) * " s.\n")
        # Storage for photon number computations.
        htPairs = Array{Float64,1}(undef, size(lPairs)[2])
        # Compute number of generated photons using heat transfer function.
        # Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
        stats = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
        write(fileStream, "Generated excitation profile computed " * string(stats.time) * " s.\n")
        write(fileStream,string(xMinN)*" "*string(xMinP)*" "*string(xMinSub)*" "*string(xMinProt)*" xMin for mbox meshing. InAs, Q, Substrate, Prot respectively. \n")
        write(fileStream,string(divNcell)*" "*string(divPcell)*" "*string(divSubCell)*" "*string(divProtCell)*" of slices for InAs, Q, and Substrate, Prot respectively.\n\n")
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




        ###### simulation for photon number
        ##simulation run with no subdivisions for accurate absolute values
        xMinProt = 0.0
        xMinN = 0.0
        xMinP = 0.0
        xMinSub = substrate_d
        divProtCell = 1
        divNcell = 1
        divPcell = 1
        divSubCell = 1
        mboxProt = 0.0
        mboxN = 0.0
        mboxP = 0.0
        mboxSub = 1.0
        #energy interval
        enrRng = (0.35, 1.0)
        # Builds slab structure, generating lVar and lPairs.
        # Arguments are: temperature of the emitter, background temperature, divisions of NCell
        stats = @timed (lVar, lPairs) = uOttawaSlabs_v3(Rad_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, firstgap_d, Rad_d, gap_d, fsf_d, emitter_d, base_d, substrate_d, Rad_dop, fsf_dop, emitter_dop, base_dop, base_As, fsf_As, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub)
        # Storage for photon number computations.
        htPairs = Array{Float64,1}(undef, size(lPairs)[2])
        # Compute number of generated photons using heat transfer function.
        # Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
        stats = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
        #copy ouput into new array to use later
        htPairsNorm = copy(htPairs)

        xMinProt = 0.0
        xMinN = 0.0
        xMinP = 0.0
        xMinSub = 0.1
        divProtCell = ceil(Int,sqrt(fsf_d*100))
        divNcell = ceil(Int,sqrt(emitter_d*100))
        divPcell = ceil(Int,sqrt(base_d*100))
        divSubCell = 30
        mboxProt = 0.0
        mboxN = 0.0
        mboxP = 0.0
        mboxSuber(x) = mboxish(x,divSubCell,substrate_d/xMinSub)
        mboxSub = find_zero(mboxSuber,2.0,Order5())
        # Open output file.
        fileStream = open("Photon Number/"*fName,"w")
        #end
        # Double check thread initialization.
        write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
        ## External program settings, internal setting contained in uOttawaSlabStruct.jl
        enrRng = (0.35, 1.0)
        ## Begin program execution.
        write(fileStream, "Execution:\n")
        # Builds slab structure, generating lVar and lPairs.
        # Arguments are: temperature of the emitter, background temperature, divisions of NCell
        # (Resolution of photon generation is thckInAs / divNCell), division of PCell
        stats = @timed (lVar, lPairs) = uOttawaSlabs_v3(Rad_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, firstgap_d, Rad_d, gap_d, fsf_d, emitter_d, base_d, substrate_d, Rad_dop, fsf_dop, emitter_dop, base_dop, base_As, fsf_As, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub)
        write(fileStream, "Layer structure constructed in " * string(stats.time) * " s.\n")
        # Storage for photon number computations.
        htPairs = Array{Float64,1}(undef, size(lPairs)[2])
        # Compute number of generated photons using heat transfer function.
        # Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
        stats = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
        write(fileStream, "Generated excitation profile computed " * string(stats.time) * " s.\n")
        write(fileStream,string(xMinN)*" "*string(xMinP)*" "*string(xMinSub)*" "*string(xMinProt)*" xMin for mbox meshing. InAs, Q, Substrate, Prot respectively. \n")
        write(fileStream,string(divNcell)*" "*string(divPcell)*" "*string(divSubCell)*" "*string(divProtCell)*" of slices for InAs, Q, and Substrate, Prot respectively.\n\n")
        # Write results to file.
        write(fileStream, "Depth (microns)	"*" Excitation Density Rate (cm-2 microns-1 s-1)\n\n")
        idxNorm = 1
        extremes = [1,divProtCell]
        # Output calculation results.
        for ind = 1 : length(htPairs)
            write(fileStream, string(round(htPairs[ind]/(sum(htPairs[extremes[1]:extremes[2]])*(lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind]-1])*evJ)*htPairsNorm[idxNorm],sigdigits=4)) * " \n")
            #to go through all the layers
            if ind == divProtCell || ind == divProtCell + divNcell || ind == divProtCell + divNcell + divPcell
                idxNorm += 1 
            end
            if ind == divProtCell
                extremes = [ind+1,ind + divNcell]
            elseif ind == divProtCell + divNcell
                extremes = [ind+1, ind + divPcell]
            elseif ind == divProtCell + divNcell + divPcell
                extremes = [ind+1,ind+divSubCell]
            end
        end
        # Flush file stream.
        close(fileStream)
    end
end




# DO NOT NEED TO CHANGE

#function to set up the 2 parameter problem as a matrixsw
@everywhere function twoParams(i::Int)
    twosweep = ceil(Int,i/length(sweep1))
    heatLayersFunc(sweep1[i - length(sweep1)*(twosweep-1)],sweep2[twosweep])
end


# perform simulations and keep track of how much there is left
@showprogress pmap(i -> twoParams(i),1:length(sweep1)*length(sweep2))



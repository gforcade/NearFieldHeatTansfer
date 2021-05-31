# - v2: to  calculate the total depth resolved heat transfer to cell

@everywhere push!(LOAD_PATH,pwd())
@everywhere begin
        
    using Distributed
    using Base.Threads, ProgressMeter, Plots, LaTeXStrings, Roots
    using FilmDataStructures, FilmUtilities, ResponseModels, OttawaSlabCell_v3_1, OttawaSlabCell_v2_1
    const evJ = 1.6021774232052327e-19
end




# (CHANGE) sweeping over 2 parameters. Change parameters in heatLayersFunc as needed
@everywhere sweep1 =  LinRange(0.01,0.2,20) 
@everywhere sweep2 = LinRange(0.4,1.0,7)#[3.0, 6.0, 10.0, 30.0, 60.0, 100.0, 300.0, 600.0, 1000.0]*(10.0^21) 


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
    Radiator_T = 700.0      #temperature of Radiator
    d_gap = 0.1             #gap between radiator and cell
    d_prot= var1 #0.05      #top cell layer thickness
    d_InAs_total = 5.0      #Emitter layer thickness
    d_InAsSbP_base = 1.5    #base layer layer thickness
    d_InAs_sub = 120.0      #substrate layer thickness
    ndoping_InAs = 3.0*(10.0^21)    #emitter layer doping
    quat_x = 0.4            #base layer As content (lattice matched to InAs)
    prot_quat_x = var2      #top layer As content (lattice matched to InAs)
    prot_ndoping = 5.0*(10.0^24)    #top layer doping (only used if layer is pure InAs)
    pdoping_InAs = 3.0*(10.0^24)    #substrate doping
    xMinProt = 0.005                #minimum top layer thickness
    xMinN = 0.05                    #minimum emitter layer thickness
    xMinP = 0.05                    #minimum base layer thickness
    xMinSub = 0.05                  #minimum substrate layer thickness
    #number of sections for a given layer
    divProtCell = ceil(Int,sqrt(d_prot*200))
    divNcell = ceil(Int,sqrt(d_InAs_total*100))
    divPcell = ceil(Int,sqrt(d_InAsSbP_base*100))
    divSubCell = 30
    #functions to get mbox which controls the section thicknesses such that they increase with depth into layer
    mboxProter(x) = mboxish(x,divProtCell,d_prot/xMinProt)
    mboxProt = find_zero(mboxProter,2.0,Order5())
    mboxNer(x) = mboxish(x,divNcell,d_InAs_total/xMinN)
    mboxN = find_zero(mboxNer,2.0,Order5())
    mboxPer(x) = mboxish(x,divPcell,d_InAsSbP_base/xMinP)
    mboxP = find_zero(mboxPer,2.0,Order5())
    mboxSuber(x) = mboxish(x,divSubCell,d_InAs_sub/xMinSub)
    mboxSub = find_zero(mboxSuber,2.0,Order5())
    
    
    
    
    # Open output file.
    #fileStream = open("photonProfile"*string(divNcell)*".txt","w")
    if prot_quat_x != 1.0
        fileStream = open("Results/Rad"*string(trunc(Int,round(Radiator_T)))*"SiGap"*string(trunc(Int,round(d_gap*1000)))*"fsf"*string(trunc(Int,round(d_prot*1000)))*"As"*string(round(prot_quat_x,sigdigits=1))*"InAs"*string(trunc(Int,round(d_InAs_total*1000)))*"Dop"*string(round(ndoping_InAs/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(d_InAsSbP_base*1000)))*"As"*string(round(quat_x,sigdigits=1))*"Sub"*string(trunc(Int,round(d_InAs_sub*1000)))*"Dop"*string(round(pdoping_InAs/1000000,sigdigits=1))*".txt","w")
    else
        fileStream = open("Results/Rad"*string(trunc(Int,round(Radiator_T)))*"SiGap"*string(trunc(Int,round(d_gap*1000)))*"fsf"*string(trunc(Int,round(d_prot*1000)))*"As"*string(round(prot_quat_x,sigdigits=1))*"Dop"*string(round(prot_ndoping/1000000,sigdigits=1))*"InAs"*string(trunc(Int,round(d_InAs_total*1000)))*"Dop"*string(round(ndoping_InAs/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(d_InAsSbP_base*1000)))*"As"*string(round(quat_x,sigdigits=1))*"Sub"*string(trunc(Int,round(d_InAs_sub*1000)))*"Dop"*string(round(pdoping_InAs/1000000,sigdigits=1))*".txt","w")
    end
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
    if prot_quat_x != 1.0
        stats = @timed (lVar, lPairs) = uOttawaSlabs_v2_1(Radiator_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, d_gap, d_prot, d_InAs_total, d_InAsSbP_base, d_InAs_sub, ndoping_InAs, pdoping_InAs, quat_x, prot_quat_x, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub)
    else
        stats = @timed (lVar, lPairs) = uOttawaSlabs_v3_1(Radiator_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, d_gap, d_prot, d_InAs_total, d_InAsSbP_base, d_InAs_sub, prot_ndoping, ndoping_InAs, pdoping_InAs, quat_x, prot_quat_x, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub)
    
    end

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
    write(fileStream, "Layer boundary edge (microns)	"*" Excitation Density Rate (W cm-2 microns-1 s-1)\n\n")
    # Location of cell boundary within the slab structure.
    #cellSurfLoc = lVar.bdrLoc[lPairs[2,1]-1]
    # Output calculation results.
    for ind = 1 : length(htPairs)-1
        #htPairs[ind] = htPairs[ind]*(lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind-1]])
        #write(fileStream, string(round(lVar.bdrLoc[lPairs[2,ind]] - cellSurfLoc, sigdigits = 4)) * "," * string(round(/(htPairs[ind], evJ), sigdigits = 4)) * "\n")
        write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(htPairs[ind]/((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind]-1])),sigdigits=4)) * " \n")
    end
    #gold backing absorption
    write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,length(htPairs)-1]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(htPairs[length(htPairs)],sigdigits=4)) * " \n")
    # Flush file stream.
    close(fileStream)
end




# DO NOT NEED TO CHANGE

#function to set up the 2 parameter problem as a matrix
@everywhere function twoParams(i::Int)
    twosweep = ceil(Int,i/length(sweep1))
    heatLayersFunc(sweep1[i - length(sweep1)*(twosweep-1)],sweep2[twosweep])
end


# perform simulations and keep track of how much there is left
@showprogress pmap(i -> twoParams(i),1:length(sweep1)*length(sweep2))

#@distributed for i = 1:length(sweep1)*length(sweep2)
#    twoParams(i)
#end

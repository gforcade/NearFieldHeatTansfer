# - v2: calculates the heat transfer vs photon energy (output is binned)

push!(LOAD_PATH,pwd())

        
using Base.Threads, ProgressMeter, LaTeXStrings, Roots
using FilmDataStructures, FilmUtilities, ResponseModels, OttawaSlabCell_v3_1, OttawaSlabCell_v2, OttawaSlabCell_v3, OttawaSlabCell_v2_1
const evJ = 1.6021774232052327e-19








function heatLayersFunc(var1::Float64,var2::Float64)
    
    #gets the mbox value for a given N=# fo  slices and r = xthck/xmin
    function mboxish(x,N,r)
        
        tot = 0.0
        for ind = 1 : N - 1
            tot += x^ind
        end
        return tot - (r - 1)
    end

    # T [K], d [um], dop [m-3]
    Radiator_T = var1 #700.0
    d_firstgap = 0.5
    d_Rad = 20.0
    d_gap = 0.1
    d_prot= 0.07
    d_InAs_total = 5.0
    d_InAsSbP_base = 1.5
    d_InAs_sub = 120.0
    Si_ndoping = 1.0*(10.0^25)
    ndoping_InAs = 3.0*(10.0^21)
    quat_x = 0.4  #smallest is best because of higher bandgap
    prot_quat_x = var2 #1.0
    prot_ndoping = 5.0*(10.0^24)
    pdoping_InAs = 3.0*(10.0^24)  #highest is best (proven)
    xMinProt = 0.0
    xMinN = 0.0
    xMinP = 0.0
    xMinSub = 0.0
    divProtCell = 1#ceil(Int,sqrt(d_prot*200))
    divNcell = 5# ceil(Int,sqrt(d_InAs_total*100))
    divPcell = 1 #ceil(Int,sqrt(d_InAsSbP_base*100))
    divSubCell = 1
    #mboxProter(x) = mboxish(x,divProtCell,d_prot/xMinProt)
    mboxProt = 0.0#find_zero(mboxProter,2.0,Order5())
    #mboxNer(x) = mboxish(x,divNcell,d_InAs_total/xMinN)
    mboxN = 0.0#find_zero(mboxNer,2.0,Order5())
    #mboxPer(x) = mboxish(x,divPcell,d_InAsSbP_base/xMinP)
    mboxP = 0.0#find_zero(mboxPer,2.0,Order5())
    #mboxSuber(x) = mboxish(x,divSubCell,d_InAs_sub/xMinSub)
    mboxSub = 0.0#find_zero(mboxSuber,2.0,Order5())
    
    ####
    ####simulations for total energy transfer
    # Open output file.
    fName = "Spectral_Resolution/Rad"*string(trunc(Int,round(Radiator_T)))*"fstGap"*string(trunc(Int,round(d_firstgap*1000)))*"Si"*string(trunc(Int,round(d_Rad*1000)))*"Dop"*string(round(Si_ndoping/1000000,sigdigits=1))*"Gap"*string(trunc(Int,round(d_gap*1000)))*"fsf"*string(trunc(Int,round(d_prot*1000)))*"As"*string(round(prot_quat_x,sigdigits=1))*"Dop"*string(round(prot_ndoping/1000000,sigdigits=1))*"InAs"*string(trunc(Int,round(d_InAs_total*1000)))*"Dop"*string(round(ndoping_InAs/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(d_InAsSbP_base*1000)))*"As"*string(round(quat_x,sigdigits=1))*"Dop"*string(round(pdoping_InAs/1000000,sigdigits=1))*"Sub"*string(trunc(Int,round(d_InAs_sub*1000)))*"Dop"*string(round(pdoping_InAs/1000000,sigdigits=1))*".txt"
    fileStream = open(fName,"w")
    # Double check thread initialization.
    write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
    write(fileStream,string(xMinN)*" "*string(xMinP)*" "*string(xMinSub)*" "*string(xMinProt)*" xMin for mbox meshing. InAs, Q, Substrate, Prot respectively. \n")
    write(fileStream,string(divNcell)*" "*string(divPcell)*" "*string(divSubCell)*" "*string(divProtCell)*" of slices for InAs, Q, and Substrate, Prot respectively.\n")
    # building simulations structure 
    (lVar, lPairs) = uOttawaSlabs_v3_1(Radiator_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, d_firstgap, d_Rad, d_gap, d_prot, d_InAs_total, d_InAsSbP_base, d_InAs_sub, Si_ndoping, prot_ndoping, ndoping_InAs, pdoping_InAs, quat_x, prot_quat_x, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub)
    ##  energy interval
    enrN = 300
    enrRngLims = LinRange(0.01, 1.0,enrN)
    
    stats = @time for i = 1:enrN-1

        enrRng = (enrRngLims[i], enrRngLims[i+1])
        # Storage for energy transfer computations.
        htPairs = Array{Float64,1}(undef, size(lPairs)[2])
        # Compute heat transfer, using heat transfer function.
        heatTfr!(lVar, lPairs, enrRng, htPairs)
        #write cell sectioning parameters in heading
        # Output calculation results.
        println(sum(htPairs))
        # Write results to file.
        write(fileStream, "\n Energy interval limits (lower and upper): "* string(enrRng[1])*" "*string(enrRng[2])*" \n")
        write(fileStream, "Layer boundary edge (microns)	"*" Excitation Density Rate (W cm-2) \n")
        for ind = 1 : length(htPairs) - 1
            write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(htPairs[ind],sigdigits=4)) * " \n")
        end
        #gold backing absorption
        write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,length(htPairs)-1]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(htPairs[length(htPairs)],sigdigits=4)) * " \n")
    end
    # Flush file stream.
    close(fileStream)
    println(stats)
    #



    #=
    ###### simulation for photon number
    xMinProt = 0.0
    xMinN = 0.0
    xMinP = 0.0
    xMinSub = 10.0
    divProtCell = 1#ceil(Int,sqrt(d_prot*200))
    divNcell = 1# ceil(Int,sqrt(d_InAs_total*100))
    divPcell = 1 #ceil(Int,sqrt(d_InAsSbP_base*100))
    divSubCell = 3
    #mboxProter(x) = mboxish(x,divProtCell,d_prot/xMinProt)
    mboxProt = 0.0#find_zero(mboxProter,2.0,Order5())
    #mboxNer(x) = mboxish(x,divNcell,d_InAs_total/xMinN)
    mboxN = 0.0#find_zero(mboxNer,2.0,Order5())
    #mboxPer(x) = mboxish(x,divPcell,d_InAsSbP_base/xMinP)
    mboxP = 0.0#find_zero(mboxPer,2.0,Order5())
    mboxSuber(x) = mboxish(x,divSubCell,d_InAs_sub/xMinSub)
    mboxSub = find_zero(mboxSuber,2.0,Order5())
    # Open output file.
    #fileStream = open("photonProfile"*string(divNcell)*".txt","w")
    fileStream = open("Results_photonNumber/Rad"*string(trunc(Int,round(Radiator_T)))*"fstGap"*string(trunc(Int,round(d_firstgap*1000)))*"Si"*string(trunc(Int,round(d_Rad*1000)))*"Dop"*string(round(Si_ndoping/1000000,sigdigits=1))*"Gap"*string(trunc(Int,round(d_gap*1000)))*"fsf"*string(trunc(Int,round(d_prot*1000)))*"As"*string(round(prot_quat_x,sigdigits=1))*"Dop"*string(round(prot_ndoping/1000000,sigdigits=1))*"InAs"*string(trunc(Int,round(d_InAs_total*1000)))*"Dop"*string(round(ndoping_InAs/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(d_InAsSbP_base*1000)))*"As"*string(round(quat_x,sigdigits=1))*"Dop"*string(round(pdoping_InAs/1000000,sigdigits=1))*"Sub"*string(trunc(Int,round(d_InAs_sub*1000)))*"Dop"*string(round(pdoping_InAs/1000000,sigdigits=1))*".txt","w")
    # Double check thread initialization.
    write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
    # Conversion of photon energy from electron volts to Joules.
    #const evJ = 1.6021774232052327e-19
    ## External program settings, internal setting contained in uOttawaSlabStruct.jl
    enrRng = (0.35, 1.3)
    ## Begin program execution.
    write(fileStream, "Execution:\n")
    # Builds slab structure, generating lVar and lPairs.
    # Arguments are: temperature of the emitter, background temperature, divisions of NCell
    # (Resolution of photon generation is thckInAs / divNCell), division of PCell
    #stats = @timed (lVar, lPairs) = uOttawaSlabs(800.0, 300.0, divNcell, 30)
    #if prot_quat_x != 1.0
    #    stats = @timed (lVar, lPairs) = uOttawaSlabs_v2(Radiator_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, d_gap, d_prot, d_InAs_total, d_InAsSbP_base, d_InAs_sub, ndoping_InAs, pdoping_InAs, quat_x, prot_quat_x, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub)
    #else
    stats = @timed (lVar, lPairs) = uOttawaSlabs_v3(Radiator_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, d_firstgap, d_Rad, d_gap, d_prot, d_InAs_total, d_InAsSbP_base, d_InAs_sub, Si_ndoping, prot_ndoping, ndoping_InAs, pdoping_InAs, quat_x, prot_quat_x, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub)
    #end
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
    # Location of cell boundary within the slab structure.
    #cellSurfLoc = lVar.bdrLoc[lPairs[2,1]-1]
    # Output calculation results.
    for ind = 1 : length(htPairs)
        #htPairs[ind] = htPairs[ind]*(lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind-1]])
        #write(fileStream, string(round(lVar.bdrLoc[lPairs[2,ind]] - cellSurfLoc, sigdigits = 4)) * "," * string(round(/(htPairs[ind], evJ), sigdigits = 4)) * "\n")
        write(fileStream, string(round(htPairs[ind]/((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind]-1])*evJ),sigdigits=4)) * " \n")
    end
    # Flush file stream.
    close(fileStream) 
    =#
end




#heatLayersFunc(700.0,0.8)
heatLayersFunc(700.0,1.0)


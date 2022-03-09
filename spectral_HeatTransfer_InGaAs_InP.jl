# - calculates the heat transfer vs photon energy (output is binned)

push!(LOAD_PATH,pwd())

        
using Base.Threads, ProgressMeter, LaTeXStrings, Roots,Dates
using FilmDataStructures, FilmUtilities, ResponseModels, OttawaSlabCell_v3_InGaAs_InP
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
    Radiator_T = 1270.0
    d_firstgap = 0.01
    d_Rad = 60.0
    d_gap = 0.1
    d_prot= 0.2
    d_InAs_total = 1.0
    d_InAsSbP_base = 0.1
    d_InAs_sub = 0.1
    Si_ndoping = -3.0*(10.0^25) 
    prot_ndoping = -1.0*(10.0^24)
    ndoping_InAs = 1.0*(10.0^23)
    pdoping_InAs = 1.0*(10.0^24) 
    substrate_dop = 1.0*(10.0^24)
    quat_x = 0.47  
    prot_quat_x = 0.47
    
    #output folder
    fName = "Results for InGaAs_InP\\"

    
    ####
    ####simulations for total energy transfer
    # Open output file.
    fName = "..\\Simulation Results\\"*fName
    fName = fName*"Spectral Heat Transfer/"*string(Dates.format(now(),"yyyy_mm_dd_H_M"))*".txt"
    fileStream = open(fName,"w")
    xMinProt = 0.0
    xMinN = 0.0
    xMinP = 0.0
    divProtCell = 1
    divNcell = 1
    divPcell = 1 
    if d_InAs_sub == 0.0 || d_InAs_sub > 30.0
        divSubCell = 0
        mboxSub = 0.0
        xMinSub = 0.0
    else
        divSubCell = 1
        mboxSub = 1.0
        xMinSub = d_InAs_sub
    end
    mboxProt = 0.0
    mboxN = 0.0
    mboxP = 0.0
    #number of layers in the simulations 
    numLayers = divNcell + divPcell + divProtCell + divSubCell

    # building simulations structure 
    (lVar, lPairs) = uOttawaSlabs_v3_InGaAs_InP(Radiator_T, 300.0, divProtCell, divNcell, divPcell, divSubCell, d_firstgap, d_Rad, d_gap, d_prot, d_InAs_total, d_InAsSbP_base, d_InAs_sub, Si_ndoping, prot_ndoping, ndoping_InAs, pdoping_InAs, substrate_dop, quat_x, prot_quat_x, mboxProt, mboxN, mboxP, mboxSub, xMinProt, xMinN, xMinP, xMinSub,0)
    # change div subcell to write it in file if thick substrate is used
    if d_InAs_sub > 30.0
        divSubCell = 1
    end
    write(fileStream,"Rad"*string(trunc(Int,round(Radiator_T)))*"fstGap"*string(trunc(Int,round(d_firstgap*1000)))*"Si"*string(trunc(Int,round(d_Rad*1000)))*"Dop"*string(round(Si_ndoping/1000000,sigdigits=1))*"Gap"*string(trunc(Int,round(d_gap*1000)))*"fsf"*string(trunc(Int,round(d_prot*1000)))*"As"*string(round(prot_quat_x,sigdigits=1))*"Dop"*string(round(prot_ndoping/1000000,sigdigits=1))*"InAs"*string(trunc(Int,round(d_InAs_total*1000)))*"Dop"*string(round(ndoping_InAs/1000000,sigdigits=1))*"Q"*string(trunc(Int,round(d_InAsSbP_base*1000)))*"As"*string(round(quat_x,sigdigits=1))*"Dop"*string(round(pdoping_InAs/1000000,sigdigits=1))*"Sub"*string(trunc(Int,round(d_InAs_sub*1000)))*"Dop"*string(round(substrate_dop/1000000,sigdigits=1))*" \n")
    write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
    write(fileStream,string(xMinN)*" "*string(xMinP)*" "*string(xMinSub)*" "*string(xMinProt)*" xMin for mbox meshing. InAs, Q, Substrate, Prot respectively. \n")
    write(fileStream,string(divProtCell)*" "*string(divNcell)*" "*string(divPcell)*" "*string(divSubCell)*" of slices for Prot, InAs, Q, and Substrate respectively.\n")
    ##  energy interval
    enrN = 300
    enrRngLims = LinRange(0.01, 1.5,enrN)

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

        for ind = 1 : numLayers
            write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(sum(htPairs[ind:numLayers+1:length(htPairs)]),sigdigits=4)) * " \n")
        end
        #gold or substrate backing absorption
        if d_InAs_sub > 30.0
            write(fileStream,  string(round(d_InAs_sub,sigdigits=4))*" "*string(round(sum(htPairs[numLayers+1:numLayers+1:length(htPairs)]),sigdigits=4)) * " \n")
            write(fileStream,  string(round(d_InAs_sub,sigdigits=4))*" "*string(round(0.0,sigdigits=4)) * " \n")
        else
            write(fileStream,  string(round((lVar.bdrLoc[lPairs[2,length(htPairs)-1]] - lVar.bdrLoc[lPairs[2,1]-1]),sigdigits=4))*" "*string(round(sum(htPairs[numLayers+1:numLayers+1:length(htPairs)]),sigdigits=4)) * " \n")
    
        end
    end
    # Flush file stream.
    close(fileStream)
    println(stats)

end




#heatLayersFunc(700.0,0.8)
heatLayersFunc(700.0,1.0)


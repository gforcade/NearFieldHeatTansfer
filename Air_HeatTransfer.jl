## heat transfer between blackbodies
## This is to validate the optical model  

#directory to code
savFileDir = pwd()#"/home/gforc034/hs_Sean_02102021_v3"

push!(LOAD_PATH,savFileDir)
using Base.Threads, ProgressMeter, Roots, Base, Printf, Plots
using FilmDataStructures, FilmUtilities, ResponseModels
const evJ = 1.6021774232052327e-19


# T [K], d [um]
Si1_T=700.0
gap_d_range = exp10.(range(log10(0.001), stop=log10(20.0), length=1000))
Si2_T=300.0

htTransfer = Float64[]

for gap_d in gap_d_range





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

    # Construct silicon response model.
    siRspE_top(enr) = cstRsp(1.0 + im * ^(10.0, -2), enr)
    siAbsE_top(enr) = imag(siRspE_top(enr))
    #Si for bottom layer
    siRspE(enr) = cstRsp(1.0 + im * ^(10.0, -2), enr)
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

    enrRng = (6.6e-4, 2.0) #energy interval
    # Builds slab structure, generating lVar and lPairs.
    # Arguments are: temperature of the emitter, background temperature, divisions of NCell



    # Storage for heat transfer computations.
    htPairs = Array{Float64,1}(undef, size(lPairs)[2])
    # Compute heat transfer using heat transfer function.
    heatTfr!(lVar, lPairs, enrRng, htPairs)

    #println("Heattransfer execution time = ",stats2.time)
    #println("Total heat transfer (W/cm2) = ", @sprintf("%.2E",sum(htPairs)))
    push!(htTransfer, sum(htPairs))
    
end

println(htTransfer)
plot(gap_d_range,htTransfer,  xaxis=:log,yaxis=:log, xlabel = "Gap (um)", ylabel = "Heat transfer (W/cm2)", thickness_scaling = 2) 









## heat transfer between doped silicon slabs
## This is to validate the optical model  

#directory to code
savFileDir = pwd()#"/home/gforc034/hs_Sean_02102021_v3"

push!(LOAD_PATH,savFileDir)
using Base.Threads, ProgressMeter, Roots, Base, Printf,CSV,DataFrames,Plots
using FilmDataStructures, FilmUtilities, ResponseModels
const evJ = 1.6021774232052327e-19


#get data to compare
data_Exp = CSV.read("../Si model validation/SiHeatTransfer_fig6_1e18.csv",DataFrame;header=2,type=Float64)
p=plot(data_Exp[!,1],data_Exp[!,2],seriestype= :scatter,label="Them",markersize=3)



# T [K], d [um], dop [cm-3]
Si1_T=1000.0
gap_d=exp10.(range(-3.0,stop=1.0,length=100))
Si2_T=300.0
Si1_dop = 1.0e10
Si2_dop = 1.0e18


    
#fileStream = open("SiHeatTransfered.txt","w")

function SiHeatTransfer(Si1_T,Si2_T,Si1_dop,Si2_dop,gap_d)
####running simulations

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
# Storage for heat transfer computations.
htPairs = Array{Float64,1}(undef, size(lPairs)[2])
# Compute heat transfer using heat transfer function.
stats2 = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)

println("Heattransfer execution time = ",stats2.time)
println("Total heat transfer (W/m2) = ", @sprintf("%.2E",sum(htPairs)*1e4))

return sum(htPairs)*1e4
end



# run function and store output
htTransfer = Array{Float64,1}(undef, size(gap_d))
for i = 1:length(gap_d)
        htTransfer[i] = SiHeatTransfer(Si1_T,Si2_T,Si1_dop,Si2_dop,gap_d[i])
end

plot!(p,gap_d,htTransfer, xlims=(1e-3,1e1), ylims=(1e4,1e8), xaxis=:log10, yaxis=:log10,     label="Us", xlabel = "Seperation distance (um)", ylabel = "Heat (W m-2)", thickness_scaling = 2, linewidth=2,legend=:topright) 
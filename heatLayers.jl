push!(LOAD_PATH,pwd())
using Base.Threads, ProgressMeter, Plots, LaTeXStrings
using FilmDataStructures, FilmUtilities, ResponseModels, OttawaSlabCell
# Open output file.
fileStream = open("photonProfile6.txt","w")
# Double check thread initialization.
write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
# Conversion of photon energy from electron volts to Joules.
const evJ = 1.6021774232052327e-19
## External program settings, internal setting contained in uOttawaSlabStruct.jl
enrRng = (0.06, 1.0)
## Begin program execution.
write(fileStream, "Execution:\n")
# Builds slab structure, generating lVar and lPairs.
# Arguments are: temperature of the emitter, background temperature, divisions of NCell
# (Resolution of photon generation is thckInAs / divNCell), division of PCell
stats = @timed (lVar, lPairs) = uOttawaSlabs(800.0, 300.0, 6, 6)
write(fileStream, "Layer structure constructed in " * string(stats.time) * " s.\n")
# Storage for photon number computations.
htPairs = Array{Float64,1}(undef, size(lPairs)[2])
# Compute number of generated photons using heat transfer function.
# Switch from flux values achieved through prefactors declared in uOttawaSlabStruct.jl.
stats = @timed heatTfr!(lVar, lPairs, enrRng, htPairs)
write(fileStream, "Generated excitation profile computed " * string(stats.time) * " s.\n\n")
# Write results to file.
write(fileStream, "Depth (microns)	"*" Excitation Density Rate (cm-2 microns-1 s-1)\n\n")
# Location of cell boundary within the slab structure.
cellSurfLoc = lVar.bdrLoc[lPairs[2,1]-1]
# Output calculation results.
for ind = 1 : length(htPairs)
	write(fileStream, string(round((lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind]-1]),digits=4)) * "    " * string(round(htPairs[ind]*(lVar.bdrLoc[lPairs[2,ind]] - lVar.bdrLoc[lPairs[2,ind]-1])/evJ,digits=4)) * "				" * string(round(/(htPairs[ind], evJ), sigdigits = 4)) * "\n")
end
# Flush file stream.
close(fileStream)

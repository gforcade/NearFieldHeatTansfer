using Base.Threads, ProgressMeter, Plots, Cubature, LaTeXStrings
using FilmDataStructures, ResponseModels, FilmUtilities, OttawaSlabCell

# Relative tolerance for cubature.
const cubRelTol = 1.0e-4

# Open output file. 
fileStream = open("dosTest.txt","w")
# Double check thread initialization. 
write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
# Conversion of photon energy from electron volts to Joules.
const evJ = 1.602e-19
## External program settings, internal setting contained in uOttawaSlabStruct.jl
enrRng = 0.1:0.01:1.5
## Begin program execution.
write(fileStream, "Execution:\n")
# Builds slab structure, generating lVar and lPairs.
# Arguments are: temperature of the emitter, background temperature, divisions of NCell 
# (Resolution is thckInAs / divNCell), division of PCell
numDivNCell = 5
numDivPCell = 5
numDivBck = 5

stats = @timed (lVar, lPairs) = uOttawaSlabs(800.0, 300.0, numDivNCell, numDivPCell, numDivBck)
write(fileStream, "Layer structure constructed in " * string(stats.time) * " s.\n")
# Storage for relative density of states computations. 
relDOSPts = Array{Float64,2}(undef, size(lPairs)[2], length(enrRng))
# Compute spatially averaged relative density of states for all layers appearing as layer pair targets. 
stats = @timed relDOSIntEval!(lVar, lPairs, enrRng, relDOSPts)
# Output results
write(fileStream, "Spectral relative density of states computed in " * string(stats.time) * " s.\n\n")
# Write results to file.
write(fileStream, "Energy (eV) "*"   Relative Transverse DOS --->\n\n")

for enrInd = 1 : length(enrRng)
	
	write(fileStream, string(round(enrRng[enrInd], sigdigits = 4)) * "		") 	

	for lyrInd = 1 : size(lPairs)[2]
	
		write(fileStream, string(round(1.0 + relDOSPts[lyrInd,enrInd], sigdigits = 4)) * "	") 
	end
	
	write(fileStream, "\n") 		
end
# Flush file stream.
close(fileStream)
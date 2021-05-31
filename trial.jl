push!(LOAD_PATH,pwd())
using Base.Threads

for i in 1:10
    println(i)
end
println(1+2)

fileStream = open("Testing.txt","w")
# Double check thread initialization.
write(fileStream, "Julia initialized with "*string(nthreads())*" threads.\n\n")
close(fileStream)
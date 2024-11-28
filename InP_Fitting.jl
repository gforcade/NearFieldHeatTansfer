# calculates the permittivities of materials and plots them

push!(LOAD_PATH,pwd())

        

using CSV,DataFrames,Plots
using eps_InAsSbP_v3,eps_InP,eps_InGaAs,eps_InAs_ntype_v2,ResponseModels, eps_gold, epsAlN




data_Exp = CSV.read("../Si model validation/n_300K_1e18.csv",DataFrame;header=0,type=Float64)


#p=plot(data_Exp[!,1],data_Exp[!,2],seriestype= :scatter,label="experimental",markersize=3)
#plot!(p,data_Exp[!,1],data_Exp[!,2])
#display(p)
#gui()



#doping in m-3


#InAs
InAs_param =eps_InAs_struct(1e20,300.0)
eps(enr) = eps_InAsntype(enr,InAs_param)

#InAsSbP
InAsSbPstructure = InAsSbP_struct(1.0,0.311*(1.0 - 1.0),1.0e21,300.0)
#eps(enr) = eps_InAsSbP_xy_ntype(enr,InAsSbPstructure)

#InGaAs
InGaAsstructure = InGaAs_struct(1.0e15,300.0)
#eps(enr) = eps_InGaAs_ntype(enr,InGaAsstructure)

#InP
InPstructure = InP_struct(1.0e23,300.0)
#eps(enr) = eps_InP_ntype(enr,InPstructure)

#AlN
AlNstructure = AlN_struct(300.0)
#eps(enr) = epsAlN_func(enr,AlNstructure)

#Si
acpDpt = dptDsc(0.045, 0.0)
dnrDpt = dptDsc(0.0456, 1e18)
siModParamsE = prmMSi(300.0, dnrDpt, acpDpt)
#eps(enr) = siRsp(enr, siModParamsE)


#gold
#eps(enr) = epsgold(enr)

#dom = Domain(data_Exp[!,1])

#photon energy
x = LinRange(0.001,1.5,1000)

#model = Model(:comp1 => FuncWrap(eps_InAsSbP_func(trial[i]),))
re = Float64[]
imre = Float64[]
@time for i = 1 : length(x)
    push!(re,real(eps(x[i])))
    push!(imre,imag(eps(x[i])))
end

# n k alpha values
n = real((re + im*imre).^0.5)
k = imag((re + im*imre).^0.5)
alph = 4.0*pi*k .* x / (1.24 * 1e-4)   #1/cm

#convert to wavelength [um]
x = 1.24 ./ x

#convert to radial frequency
#x =  x / 6.582e-16 



#=
for i = 1 : length(LinRange(0.01,1,100))
    push!(re,real(eps_InAsSbP_func(trial[i])))
    push!(imer,imag(eps_InAsSbP_func(trial[i])))
end
=#

#plot!(p,x,n, xlims=(1,1e2), ylims=(1,5), xaxis=:log10,      label="model", xlabel = "Wavelength (um)", ylabel = "Refractive Index", thickness_scaling = 2, linewidth=2,legend=:bottomleft) #xlims = (0.0,2.0),ylims = (0.1,10.0)
#display(p) 
#sleep(10)
#println(re)
#println(imer)




#write nk data to file
fileStream = open("InAs_6e14.txt","w")
write(fileStream,  "*Wavelength (um) n k \n")
write(fileStream,"ComplexRefractiveIndex { \n")
write(fileStream,"  Formula = 1 \n")
write(fileStream, "NumericalTable( \n")

x = reverse(x)
n=reverse(n)
k=reverse(k)

for i = 1 : length(x)

    write(fileStream,  string(round(x[i],sigdigits=4))*" "*string(round(n[i],sigdigits=4))*" "*string(round(k[i],sigdigits=4)) * "; \n")
end
write(fileStream, ") \n  } \n")



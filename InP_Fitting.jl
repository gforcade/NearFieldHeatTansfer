# calculates the permittivities of materials and plots them

push!(LOAD_PATH,pwd())

        

using CSV,DataFrames,Plots
using eps_InAsSbP_v3,eps_InP,eps_InGaAs,eps_InAs_ntype_v2,ResponseModels, eps_gold




data_Exp = CSV.read("Eps1_exp_InAs.csv",DataFrame;header=2,type=Float64)


#p=plot(data_Exp[!,1],data_Exp[!,2],seriestype= :scatter,label="experimental")
#plot!(p,data_Exp[!,1],data_Exp[!,2])
#display(p)
#gui()



#doping in m-3

#InAsSbP
InAsSbPstructure = InAsSbP_struct(1.0,0.311*(1.0 - 1.0),1.0e21,300.0)
#eps(enr) = eps_InAsSbP_xy_ntype(enr,InAsSbPstructure)

#InGaAs
InGaAsstructure = InGaAs_struct(1.0e23,300.0)
#eps(enr) = eps_InGaAs_ntype(enr,InGaAsstructure)

#InP
InPstructure = InP_struct(1.0e23,300.0)
#eps(enr) = eps_InP_ntype(enr,InPstructure)

#Si
acpDpt = dptDsc(0.045, 1e19)
dnrDpt = dptDsc(0.0456, 0.0)
siModParamsE = prmMSi(600.0, dnrDpt, acpDpt)
eps(enr) = siRsp(enr, siModParamsE)



#eps(enr) = epsgold(enr)

#dom = Domain(data_Exp[!,1])

x = LinRange(0.001,1.2,1000)

#model = Model(:comp1 => FuncWrap(eps_InAsSbP_func(trial[i]),))
re = Float64[]
imre = Float64[]
@time for i = 1 : length(x)
    push!(re,real(eps(x[i])))
    push!(imre,imag(eps(x[i])))
end

# n k values
n = real((re + im*imre).^0.5)
k = imag((re + im*imre).^0.5)

#convert to wavelength
x = 1.24 ./ x

#convert to radial frequency
#x =  x / 6.582e-16 



#=
for i = 1 : length(LinRange(0.01,1,100))
    push!(re,real(eps_InAsSbP_func(trial[i])))
    push!(imer,imag(eps_InAsSbP_func(trial[i])))
end
=#

p =plot(x,k,xlims = (1, 100), ylims=(1e-4,10.0), xaxis=:log,yaxis=:log,  label="model", xlabel = "Wavelength (um)", ylabel = "Refractive Index", thickness_scaling = 2) #xlims = (0.0,2.0),ylims = (0.1,10.0)
display(p) 
#sleep(10)
#println(re)
#println(imer)

#==#


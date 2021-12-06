# calculates the permittivities of materials and plots them

push!(LOAD_PATH,pwd())

        

using CSV,DataFrames,Plots
using eps_InAsSbP_v3#,eps_InAs_ntype,ResponseModels




data_Exp = CSV.read("Eps1_exp_InAs.csv",DataFrame;header=2,type=Float64)


p=plot(data_Exp[!,1],data_Exp[!,2],seriestype= :scatter,label="experimental")
#plot!(p,data_Exp[!,1],data_Exp[!,2])
display(p)
gui()
sleep(10)

#

InAsSbPstructure = InAsSbP_struct(1.0,0.311*(1.0 - 1.0),1.0e21,300.0)
eps_InAsSbP_func(enr) = eps_InAsSbP_xy_ntype(enr,InAsSbPstructure)


#dom = Domain(data_Exp[!,1])

x = LinRange(0.01,1.0,100)

#model = Model(:comp1 => FuncWrap(eps_InAsSbP_func(trial[i]),))
re = []
imer = []
for i = 1 : length(x)
    push!(re,real(eps_InAsSbP_func(x[i])))
    push!(imer,imag(eps_InAsSbP_func(x[i])))
end
#=
for i = 1 : length(LinRange(0.01,1,100))
    push!(re,real(eps_InAsSbP_func(trial[i])))
    push!(imer,imag(eps_InAsSbP_func(trial[i])))
end
=#
plot!(p,x,re,xlims = (0.0,1.0), ylims = (11.0,14.0), label="model", xlabel = "Photon Energy (eV)", ylabel = "Epsilon 1", thickness_scaling = 2)
display(p)
gui()
sleep(10)
#println(re)
#println(imer)

#==#


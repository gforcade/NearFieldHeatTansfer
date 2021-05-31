# - v3: calculates the total depth resolved heat transfer + photon number to cell

push!(LOAD_PATH,pwd())

        

using CSV,DataFrames,Plots
#using eps_InAsSbP_v2#,eps_InAs_ntype,ResponseModels




data_Exp = CSV.read("Eps1_exp_InP.csv",DataFrame;header=2,type=Float64)


p=plot(data_Exp[!,1],data_Exp[!,2],seriestype= :scatter)
plot!(p,data_Exp[!,1],data_Exp[!,2])
display(p)
gui()
sleep(10)
#=
InAsSbPstructure = InAsSbP_struct(0.4,0.311*(1 - 0.4),1.0e24,300.0)
eps_InAsSbP_func(enr) = eps_InAsSbP_xy_ptype(enr,InAsSbPstructure)


dom = Domain(data_Exp[!,3])

model = Model(:comp1 => FuncWrap(eps_InAsSbP_func(trial[i]),))
re = []
imer = []
for i = 1 : length(LinRange(0.01,1,100))
    push!(re,real(eps_InAsSbP_func(trial[i])))
    push!(imer,imag(eps_InAsSbP_func(trial[i])))
end
#=
for i = 1 : length(LinRange(0.01,1,100))
    push!(re,real(eps_InAsSbP_func(trial[i])))
    push!(imer,imag(eps_InAsSbP_func(trial[i])))
end
=#
p=plot(trial,imer)
display(p)
gui()
sleep(10)
#println(re)
#println(imer)

=#


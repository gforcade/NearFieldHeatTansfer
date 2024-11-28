# - v3: calculates the total depth resolved heat transfer + photon number to cell

push!(LOAD_PATH,pwd())

        

using Plots,DataFrames,LaTeXStrings
using eps_InAsSbP_v3,ResponseModels
const evJ = 1.6021774232052327e-19
const hbar = 1.05457*^(10,-34) 
default(show=true)

#experimental data
#using CSV
#data_Exp = CSV.read("Absorption_Coefficient_InAs.csv",DataFrame;header=2,type=Float64)#"Eps2_exp_InAs.csv"
#p=plot(data_Exp[!,1],data_Exp[!,2],seriestype= :scatter,xaxis=:log,yaxis=:log,xtickfontsize=16,ytickfontsize=16,legendfontsize=16,legend=:bottomright,framestyle=:box,label="\$\\mathrm{Experimental}\$")

@inline function eVtoOmega(energy) #energy in eV
    return(energy*evJ/hbar)
end
precompile(eVtoOmega, (Float64,))

x = 1.0
InAsSbPstructure = InAsSbP_struct(x,0.311*(1 - x),5.0e24,300.0)
#eps_InAsSbP_func(enr) = eps_InAsSbP_xy_ntype(enr,InAsSbPstructure)


dops = [8.3e18*1e6]

for i = 1:length(dops)

    InAs_param =eps_InAs_struct(dops[i],300.0)
    eps_InAsSbP_func(enr) = eps_InAsntype(enr,InAs_param)


    trial = LinRange(0.01,0.68,1000)



    sim = []
    for i = 1 : length(trial)
        eps = eps_InAsSbP_func(trial[i])
        alph = eVtoOmega(trial[i])/3e10*sqrt(2.0)*sqrt(sqrt(eps.re^2 + eps.im^2)-eps.re)
        push!(sim,alph)
    end


    #println(re)
    display(plot(trial,sim,yaxis=:log,linewidth=3,label="\$\\mathrm{n="*string(dops[i]/1e6)*"}\\mathrm{(cm}^{-3}\\mathrm{)}\$")) 
    #plot!(data_Exp[1:18,1],data_Exp[1:18,2],seriestype= :scatter)#,yaxis=:log)
end
ylims!(10,1000)
ylabel!("\$\\alpha \\; \\mathrm{(cm}^{-1}\\mathrm{)}\$",labelfontsize=16)
xlabel!("\$\\mathrm{Energy \\; (eV)}\$",labelfontsize=16)
savefig("Asborption_Coefficient.svg")

gui()

#println(re)





# - v3: calculates the total depth resolved heat transfer + photon number to cell

push!(LOAD_PATH,pwd())

        

using Plots,CSV,DataFrames,LaTeXStrings
using eps_InAsSbP_v3,ResponseModels
const evJ = 1.6021774232052327e-19
default(show=true)

#experimental data
data_Exp = CSV.read("Absorption_Coefficient_InAs.csv",DataFrame;header=2,type=Float64)#"Eps2_exp_InAs.csv"
p=plot(data_Exp[!,1],data_Exp[!,2],seriestype= :scatter,xaxis=:log,yaxis=:log,xtickfontsize=16,ytickfontsize=16,legendfontsize=16,legend=:bottomright,framestyle=:box,label="\$\\mathrm{Experimental}\$")



x = 0.4
InAsSbPstructure = InAsSbP_struct(x,0.311*(1 - x),5.0e24,300.0)
#eps_InAsSbP_func(enr) = eps_InAsSbP_xy_ntype(enr,InAsSbPstructure)


dops = [3.0e21]

for i = 1:length(dops)

    InAs_param =eps_InAs_struct(dops[i],300.0)
    eps_InAsSbP_func(enr) = eps_InAsntype(enr,InAs_param)


    trial = LinRange(0.01,0.9,1000)



    sim = []
    for i = 1 : length(trial)
        push!(sim,imag(eps_InAsSbP_func(trial[i])))
    end


    #println(re)
    plot!(trial,sim,xaxis=:log,yaxis=:log,linewidth=3,label="\$\\mathrm{n="*string(dops[i])*"}\\mathrm{(cm}^{-3}\\mathrm{)}\$")
    #plot!(data_Exp[1:18,1],data_Exp[1:18,2],seriestype= :scatter)#,yaxis=:log)
end

ylabel!("\$\\alpha \\; \\mathrm{(cm}^{-1}\\mathrm{)}\$",labelfontsize=16)
xlabel!("\$\\mathrm{Energy \\; (eV)}\$",labelfontsize=16)
#savefig("Asborption_Coefficient.svg")
display(p)
#gui()
sleep(10)
#println(re)
#println(imer)




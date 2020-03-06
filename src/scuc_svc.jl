cd("/home/cxlab/syh/repo/SRUC")
clearconsole()
using Revise
using JuMP,Gurobi,Suppressor,PowerModels,CSV,Plots
Base.GC.enable(false)
Base.GC.enable(false)
includet("Problems.jl")
wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)
# define constants and input data
const T = 24
silence()
# case = PowerModels.parse_file(joinpath(@__DIR__, "../res/1-1-2019-98-1.m"))
case = parse_file("input/case5.m")
pm = build_model(case, ACPPowerModel,PowerModels.post_opf)
case_dict = pm.ref[:nw][0];
pmax = maximum([case_dict[:gen][gen]["pmax"] for gen in keys(case_dict[:gen])]);
battery = Dict("pmax_in"=>pmax/6,"pmax_out"=>pmax/6,"η_in"=>0.7,"η_out"=>0.7,"C_max"=>pmax/6 * 6,"SOC_max"=>0.9,"SOC_min"=>0.1,"SOC_int"=>0.5,"deg_cost"=>0.01);
case_dict[:battery] = Dict([(i,battery) for i in range(1,stop=4)]);
case_dict[:windfarm] = Dict([(i,0) for i in range(1,stop=2)]);
for gen in keys(case_dict[:gen])
    # if case_dict[:gen][gen]["pmin"] == 0
    #     case_dict[:gen][gen]["pmin"] = 0.1 * case_dict[:gen][gen]["pmax"]
    # end
    if case_dict[:gen][gen]["pmax"] == 0
        delete!(case_dict[:gen],gen);
    else
        case_dict[:gen][gen]["pmin"] = case_dict[:gen][gen]["pmax"]/5;
        case_dict[:gen][gen]["startup"] = 500;
    end
end
wind_power = CSV.read("input/wind_2020.csv")[1:T,:]
load = CSV.read("input/load data.csv")[1:T,:load]
loadmax = 1.2*sum(max(case_dict[:gen][gen]["pmax"],0) for gen in keys(case_dict[:gen]))
load_real = load * loadmax/maximum(load)
r = 0.3 * loadmax/maximum(wind_power[:wf1])
data = [Dict(:wind_power=>Dict(1=>r*wind_power[t,:wf1],2=>r*wind_power[t,:wf2]),
            :load=>load_real[t]) for t in 1:T]

# problem construction
clearconsole()
dayahead = Problems.dayaheadProblem(case_dict,data)
@info("calculating initial solution by prior list...")
dayahead = Problems.prior_list_modification(dayahead,case_dict)
@info("optimizing dayahead decision...")
optimize!(dayahead)
@info("loadCut: $(value(sum(dayahead[:loadCut])))")
@info("ug: $(sum(value.(dayahead[:ug]).data,dims=1))")
# initialization
intraday = @suppress [Problems.intradayProblem(case_dict,data[t]) for t in 1:T]
intraday[T] = Problems.fix_tail(intraday[T],case_dict)# tackle the tail effect
intradayMax = deepcopy(intraday);
uncertainties = [Problems.PolygonUncertaintySet([0,0,0],[1 0 0;0 1 0;0 0 1],1.5),Problems.PolygonUncertaintySet([1,1,1],[1 0 0;0 1 0;0 0 1],1.5)]
# uncertainties = [Problems.PolygonUncertaintySet([0,0,0],[1 0 0;0 1 0;0 0 1],1.5)]
vertice = Problems.getUnionVertice(uncertainties)
vertice_series = [vertice for t in 1:T]
# @suppress_err begin
# for t in 1:T
#     optimize!(intradayMax[t],Gurobi.Optimizer(MIPGap=1e-3,OutputFlag=1))
# end
@info("intraday dispatch started\n")
intraday,intradayMax,additional= Problems.RDDP(intraday,intradayMax,dayahead,vertice_series,case_dict,data)
Ses_curve_dayahead = [[value(dayahead[:Ses][n,t]) for t in 1:T] for n in 1:4]
Ses_curve = [[value(intraday[t][:Ses][n]) for t in 1:T] for n in 1:4]
# wind_worst = [sum(value(intradayMax[t][:Pw_max][n]) for n in 1:2) for t in 1:T]
plot(Ses_curve)
real = [sum(value(intraday[τ][:cost_now]) for τ in t:T) for t in 1:T]
upper = additional[:upper]
# upper = [objective_value(intradayMax[t]) for t in 1:T]
lower = [objective_value(intraday[t]) for t in 1:T]
# plot(wind_worst)
# mean_cost,worst_cost = Problems.evalutaion(dayahead,intraday,intradayMax,case_dict,data)
# printstyled("------------SUMMARY-------------\n";color=:green)
# printstyled("Number of iteration: $n_iter\n";color=:blue)
# printstyled("Intraday expected cost: $mean_cost\n";color=:blue)
# printstyled("Intraday worst cost: $worst_cost\n";color=:blue)
Base.GC.enable(true)
Base.GC.enable(true)
Base.GC.gc()
plot([upper,lower,real])

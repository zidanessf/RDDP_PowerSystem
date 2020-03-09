cd("/home/cxlab/syh/repo/SRUC")
ENV["GKSwstype"]="100"
#clearconsole()
using Distributed,JSON
config = Dict()
open("src/CONFIG.json", "r") do f
    global config
    config=JSON.parse(f)  # parse and transform data
end
if nprocs() < config["CoreNumber"]
    addprocs(config["CoreNumber"] - nprocs())
end
@everywhere using Revise,Suppressor
using JuMP,Gurobi,PowerModels,CSV,Plots
# Base.GC.enable(false)
# Base.GC.enable(false)
if nprocs == 1
    includet("src/Problems.jl")
else
    @everywhere includet("src/Problems.jl")
end
# define constants and input data
const T = config["T"]
silence()
# case = PowerModels.parse_file(joinpath(@__DIR__, "../res/1-1-2019-98-1.m"))
casename = config["casename"]
case = parse_file("input/"*casename*".m")
pm = build_model(case, ACPPowerModel,PowerModels.post_opf)
case_dict = pm.ref[:nw][0];
pmax = maximum([case_dict[:gen][gen]["pmax"] for gen in keys(case_dict[:gen])]);
battery = config["battery"]
case_dict[:battery] = Dict([(i,battery) for i in range(1,stop=2)]);
case_dict[:windfarm] = Dict([(1,1),(2,3)]);
case_dict[:battery_location] = Dict([(1,1),(2,3)]);
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
loadmax = sum(max(case_dict[:gen][gen]["pmax"],0) for gen in keys(case_dict[:gen]))
load_real = load * loadmax/maximum(load)
r = 0.15 * loadmax/maximum(wind_power[:wf1])
data = [Dict(:wind_power=>Dict(1=>r*wind_power[t,:wf1],2=>r*wind_power[t,:wf2]),
            :load=>load_real[t]) for t in 1:T]

# problem construction
#clearconsole()
dayahead = Problems.dayaheadProblem(case_dict,data)
@info("calculating initial solution by prior list...")
dayahead = Problems.prior_list_modification(dayahead,case_dict)
@info("optimizing dayahead decision...")
optimize!(dayahead)
@info("loadCut: $(value(sum(dayahead[:loadCut])))")
@info("ug: $(sum(value.(dayahead[:ug]).data,dims=1))")
# initialization
uncertainties = [Problems.PolygonUncertaintySet([0,0],[1 0;0 1],1.5),Problems.PolygonUncertaintySet([1,1],[1 0;0 1],1.5)]
# uncertainties = [Problems.PolygonUncertaintySet([0,0,0],[1 0 0;0 1 0;0 0 1],1.5)]
vertice = Problems.getUnionVertice(uncertainties)
vertice_series = [vertice for t in 1:T]
if nprocs == 1
    RTD = @suppress Problems.RealTimeDispatchModel(
        (Dict(:ug=>value.(dayahead[:ug]),:gens=>dayahead[:gens])),
        [Problems.intradayProblem(case_dict,data[t],config["has_pf"]) for t in 1:T],
        [Problems.intradayProblem(case_dict,data[t],config["has_pf"])  for t in 1:T],
        [Dict() for t in 1:T],
        vertice_series,
        data,
        case_dict)
    Problems.fix_tail()
else
    @everywhere RTD = @suppress Problems.RealTimeDispatchModel(
        $(Dict(:ug=>value.(dayahead[:ug]),:gens=>dayahead[:gens])),
        [Problems.intradayProblem($case_dict,$data[t],config["has_pf"]) for t in 1:$T],
        [Problems.intradayProblem($case_dict,$data[t],config["has_pf"])  for t in 1:$T],
        [Dict() for t in 1:$T],
        $vertice_series,
        $data,
        $case_dict)
    @everywhere Problems.fix_tail()
end
# @suppress_err begin
# for t in 1:T
#     optimize!(intradayMax[t],Gurobi.Optimizer(MIPGap=1e-3,OutputFlag=1))
# end
@info("intraday dispatch started\n")
solution_status= Problems.RDDP(RTD)
c1,c2,c3 = Problems.evalutaion()
CSV.write("results/"*casename*"_solver_report.csv",solution_status) 
Ses_curve = [[value(Main.RTD.intraday[t][:Ses][n]) for t in 2:T] for n in keys(case_dict[:battery])]
# wind_worst = [sum(value(intradayMax[t][:Pw_max][n]) for n in 1:2) for t in 1:T]
plot(Ses_curve)
s1 = [value(RTD.intraday[t][:cost_to_go]) for t in 2:T]
s2 = [value(RTD.intradayMax[t][:cost_to_go]) for t in 2:T]
plot([s1,s2])
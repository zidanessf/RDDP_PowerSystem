cd("/home/cxlab/syh/repo/SRUC")
ENV["GKSwstype"]="100"
#clearconsole()
using Distributed
if nprocs() < 6
    addprocs(6 - nprocs())
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
const T = 24
silence()
# case = PowerModels.parse_file(joinpath(@__DIR__, "../res/1-1-2019-98-1.m"))
casename = "case5"
case = parse_file("input/"*casename*".m")
pm = build_model(case, ACPPowerModel,PowerModels.post_opf)
case_dict = pm.ref[:nw][0];
pmax = maximum([case_dict[:gen][gen]["pmax"] for gen in keys(case_dict[:gen])]);
battery = Dict("pmax_in"=>pmax/6,"pmax_out"=>pmax/6,"η_in"=>0.7,"η_out"=>0.7,"C_max"=>pmax/6 * 6,"SOC_max"=>0.9,"SOC_min"=>0.1,"SOC_int"=>0.5,"deg_cost"=>0.01);
case_dict[:battery] = Dict([(i,battery) for i in range(1,stop=3)]);
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
#clearconsole()
dayahead = Problems.dayaheadProblem(case_dict,data)
@info("calculating initial solution by prior list...")
dayahead = Problems.prior_list_modification(dayahead,case_dict)
@info("optimizing dayahead decision...")
optimize!(dayahead)
@info("loadCut: $(value(sum(dayahead[:loadCut])))")
@info("ug: $(sum(value.(dayahead[:ug]).data,dims=1))")
# initialization
uncertainties = [Problems.PolygonUncertaintySet([0,0,0],[1 0 0;0 1 0;0 0 1],1.5),Problems.PolygonUncertaintySet([1,1,1],[1 0 0;0 1 0;0 0 1],1.5)]
# uncertainties = [Problems.PolygonUncertaintySet([0,0,0],[1 0 0;0 1 0;0 0 1],1.5)]
vertice = Problems.getUnionVertice(uncertainties)
vertice_series = [vertice for t in 1:T]
if nprocs == 1
    RTD = @suppress Problems.RealTimeDispatchModel(
        (Dict(:ug=>value.(dayahead[:ug]),:gens=>dayahead[:gens])),
        [Problems.intradayProblem(case_dict,data[t]) for t in 1:T],
        [Problems.intradayProblem(case_dict,data[t])  for t in 1:T],
        [Dict() for t in 1:T],
        vertice_series,
        data,
        case_dict)
    Problems.fix_tail()
else
    @everywhere RTD = @suppress Problems.RealTimeDispatchModel(
        $(Dict(:ug=>value.(dayahead[:ug]),:gens=>dayahead[:gens])),
        [Problems.intradayProblem($case_dict,$data[t]) for t in 1:$T],
        [Problems.intradayProblem($case_dict,$data[t])  for t in 1:$T],
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
CSV.write("results/"*casename*"_solver_report.csv",solution_status) 
Ses_curve = [[value(Main.RTD.intraday[t][:Ses][n]) for t in 1:T] for n in keys(case_dict[:battery])]
# wind_worst = [sum(value(intradayMax[t][:Pw_max][n]) for n in 1:2) for t in 1:T]
plot(Ses_curve)

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
using JuMP,Gurobi,PowerModels,CSV,Plots,Random,Interpolations,DataFrames,Statistics,StatsPlots,Distributions
# Base.GC.enable(false)
# Base.GC.enable(false)
if nprocs() == 1
    includet("Problems.jl")
else
    @everywhere includet("src/Problems.jl")
end
# define constants and input data
const T = config["T"]
silence()
# case = PowerModels.parse_file(joinpath(@__DIR__, "../res/1-1-2019-98-1.m"))
casename = config["casename"]
case = PowerModels.parse_file("input/"*casename*".m")
pm = build_model(case, ACPPowerModel,PowerModels.post_opf)
case_dict = pm.ref[:nw][0];
pmax = maximum([case_dict[:gen][gen]["pmax"] for gen in keys(case_dict[:gen])]);
case_dict[:battery] = Dict([(i,config["battery"]) for i in range(1,stop=1)]);
# case_dict[:windfarm] = Dict([(1,3),(2,5)]);
# case_dict[:battery_location] = Dict([(1,3),(2,5)]);
case_dict[:windfarm] = Dict([(1,3)]);
case_dict[:battery_location] = Dict([(1,3)]);
Random.seed!(1234)
for gen in keys(case_dict[:gen])
    # if case_dict[:gen][gen]["pmin"] == 0
    #     case_dict[:gen][gen]["pmin"] = 0.1 * case_dict[:gen][gen]["pmax"]
    # end
    # case_dict[:gen][gen]["cost"][end-1] = case_dict[:gen][gen]["cost"][end-1] + 100*rand()
    if case_dict[:gen][gen]["pmax"] == 0
        delete!(case_dict[:gen],gen);
    else
        case_dict[:gen][gen]["pmin"] = case_dict[:gen][gen]["pmax"]/5;
        case_dict[:gen][gen]["startup"] = 500;
    end
end
wind_power = CSV.read("input/wind_2020.csv")[1:24,:]
load = CSV.read("input/load data.csv")[1:24,:load]
loadmax = sum(max(case_dict[:gen][gen]["pmax"],0) for gen in keys(case_dict[:gen]))
load_real = 0.8*load * loadmax/maximum(load)
r = 0.15 * loadmax/maximum(wind_power[:wf1])
load_ipt = extrapolate(interpolate(load_real,BSpline(Linear())),Flat())
wind_ipt1 = extrapolate(interpolate(wind_power[:,:wf1],BSpline(Linear())),Flat())
wind_ipt2 = extrapolate(interpolate(wind_power[:,:wf2],BSpline(Linear())),Flat())
data = [Dict(:wind_power=>Dict(1=>r*wind_ipt1(t*24/T),2=>r*wind_ipt2(t*24/T)),
            :load=>load_ipt(t*24/T)) for t in 1:T]

# problem construction
Problems.makePTDF(case_dict)

#clearconsole()
dayahead = Problems.dayaheadProblem(case_dict,data,config)
# @info("calculating initial solution by prior list...")
# Problems.prior_list_modification(dayahead,case_dict)
# @info("optimizing dayahead decision...")
# optimize!(dayahead)
# @info("loadCut: $(value(sum(dayahead[:loadCut])))")
# @info("ug: $(sum(value.(dayahead[:ug]).data,dims=1))")
# initialization
# uncertainties = [Problems.PolygonUncertaintySet([0,0],[1 0;0 1],1.5),Problems.PolygonUncertaintySet([1,1],[1 0;0 1],1.5)]
# uncertainties = [Problems.PolygonUncertaintySet([0,0,0],[1 0 0;0 1 0;0 0 1],1.5)]
uncertainties = [Problems.PolygonUncertaintySet(0,1,1)]
vertice = Problems.getUnionVertice(uncertainties)
# unit_commitment = value.(dayahead[:ug])
unit_commitment = [1 for i in dayahead[:gens]]
vertice_series = [vertice for t in 1:T]
@time if nprocs() == 1
    RTD = @suppress_out Problems.RealTimeDispatchModel(
        (Dict(:ug=>unit_commitment,:gens=>dayahead[:gens])),
        [Problems.intradayProblem(case_dict,data[t],config) for t in 1:T],
        [Problems.intradayProblem(case_dict,data[t],config)  for t in 1:T],
        [Dict() for t in 1:T],
        vertice_series,
        data,
        case_dict)
    Problems.fix_tail()
else
    @everywhere RTD = @suppress Problems.RealTimeDispatchModel(
        $(Dict(:ug=>unit_commitment,:gens=>dayahead[:gens])),
        [Problems.intradayProblem($case_dict,$data[t],$config) for t in 1:$T],
        [Problems.intradayProblem($case_dict,$data[t],$config)  for t in 1:$T],
        [Dict() for t in 1:$T],
        $vertice_series,
        $data,
        $case_dict)
    @everywhere Problems.fix_tail()
end
println("")
@info("intraday dispatch started\n")
solution_status= Problems.RDDP(RTD,config)
maxima = []
objectives,worst = Problems.evalutaion(RTD,500)
dist,μ,σ = Problems.peaksOverThresholdEstimator(objectives,0.95)
s1 = [value(RTD.intraday[t][:cost_now]) for t in 1:T]
s2 = [value(RTD.intradayMax[t][:cost_now]) for t in 1:T]
# plot(dist)
CSV.write("results/"*casename*"_solver_report.csv",solution_status) 
d = Dict()
for gen in RTD.dayahead[:gens]
    d["Unit_$gen"] = [value(RTD.intraday[t][:Pg][gen]) for t in 1:T]
end
for gen in keys(RTD.case_dict[:battery])
    d["ESS_$gen"] = [value(RTD.intraday[t][:Ses][gen]) for t in 1:T]
end
CSV.write("results/"*casename*"_decisions.csv",DataFrame(d)) 
Ses_curve = [[value(Main.RTD.intraday[t][:Ses][n]) for t in 1:T] for n in keys(case_dict[:battery])]
# wind_worst = [sum(value(intradayMax[t][:Pw_max][n]) for n in 1:2) for t in 1:T]
plot(Ses_curve)
# s1 = [value(RTD.intraday[t][:cost_to_go]) for t in 1:T]
# s2 = [value(RTD.intradayMax[t][:cost_to_go]) for t in 1:T]

# plot([s1,s2])
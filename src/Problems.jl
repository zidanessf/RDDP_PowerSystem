module Problems
export dayaheadProblem,intradayProblem
using JuMP,Gurobi,Suppressor,Random,CDDLib,Polyhedra,DataFrames,Distributed
using SharedArrays,JSON,LinearAlgebra,Statistics,Distributions,NLsolve
# config = Dict()
# open("src/CONFIG.json", "r") do f
#     global config
#     config=JSON.parse(f)  # parse and transform data
# end
# const SHARED_GUROBI_ENV = [Gurobi.Env() for i = 1:25]
@everywhere using JuMP
@everywhere using Distributed
include("print_iteration.jl")
struct RealTimeDispatchModel
    dayahead::Dict{}
    intraday::Array{JuMP.Model}
    intradayMax::Array{JuMP.Model}
    states::Array{Dict}
    vertice::Array{}
    data::Array{}
    case_dict::Dict
end
total_solves = 0
number_of_reject = 0
sceanio_not_change_num = 0
# function intradayProblemBL(ref,data,sense,cost_to_go_bound=1e8)
#     # T,ramp,UT,DT,γ_load,γ_wind,α = 24,0.25,4,4,10000,5000,0.2
#     blm = BilevelModel()
#     wst = Upper(blm)
#     m = Lower(blm)
#     aggr = []
#     gens_param_lost = [gen for gen in keys(ref[:gen]) if length(ref[:gen][gen]["cost"])<2]
#     gens = setdiff([gen for gen in keys(ref[:gen])],aggr)
#     gens = setdiff(gens,gens_param_lost)
#     ess = [gen for gen in keys(ref[:battery])]
#     wfs = [gen for gen in keys(ref[:windfarm])]
#     @variable(wst,Pw_max[wfs])
#     @variable(wst,Pw_err[wfs])
#     @variable(m,Pes[ess])
#     @variable(m,t[ess])
#     @variable(m,Ses[ess])
#     @variable(m,Sesa[ess])#ancestor
#     @variable(m,SocGuide[ess])
#     @variable(m,Pg[gens])
#     @variable(m,Pgmax[gens])
#     @variable(m,Pgu[gens])
#     @variable(m,Pgl[gens])
#     @variable(m,Pgcost[gens])
#     @variable(m,Pw[wfs])
#     @variable(m,cost_to_go)
#     @variable(m,loadCut)
#     @variable(m,unbalance)
#     #network
#     @variable(m,θ[keys(ref[:bus])])
#     #objective value
#     fuel_cost = @expression(m,sum(Pgcost))
#     deg_cost = @expression(m,sum(t[gen]*ref[:battery][gen]["deg_cost"] for gen in ess))
#     wind_Cur_cost = @expression(m,γ_wind*(sum(Pw_max) - sum(Pw)))
#     load_Cut_cost = @expression(m,γ_load*unbalance)
#     # ramp_cost = @expression(m,γ_load*sum(unable_to_ramp_down))
#     # obj = fuel_cost + deg_cost + wind_Cur_cost + load_Cut_cost
#     obj = fuel_cost + load_Cut_cost + deg_cost + wind_Cur_cost + cost_to_go + sum(SocGuide)
#     @objective(m,Min,obj)
#     if sense == "max"
#         @objective(wst,Max,obj)
#     else
#         @objective(wst,Min,obj)
#     end
#     @constraint(m,cost_to_go >= 0)
#     @constraint(m,cost_to_go <= cost_to_go_bound)
#     for gen in wfs
#         @constraint(wst,Pw_err[gen] >= -α*data[:wind_power][gen])
#         @constraint(wst,Pw_err[gen] <= α*data[:wind_power][gen])
#         @constraint(wst,Pw_max[gen] == Pw_err[gen] + data[:wind_power][gen])
#     end
#     for gen in gens
#         # @constraint(m,Pg[gen] <= ref[:gen][gen]["pmax"])#to avoid unboundness
#         # @constraint(m,Pg[gen] >= 0) #to avoid unboundness 有两个界容易算不出来，不知为何
#         @constraint(m,Pg[gen] <= Pgu[gen])
#         @constraint(m,Pg[gen] >= Pgl[gen])
#         @constraint(m,Pgcost[gen] >= 0)
#         pcostmax = ref[:gen][gen]["cost"][2]*ref[:gen][gen]["pmax"] + (ref[:gen][gen]["cost"][1])*ref[:gen][gen]["pmax"]^2
#         @constraint(m,Pgcost[gen] <= 2*pcostmax)
#         for k = 1:K_seg
#             pk = ref[:gen][gen]["pmax"]*k/K_seg
#             costk = (ref[:gen][gen]["cost"][2])*pk + (ref[:gen][gen]["cost"][1])*pk^2
#             ratek = 2*ref[:gen][gen]["cost"][1]*pk
#             @constraint(m,Pgcost[gen] >= costk + ratek*(Pg[gen] - pk))
#         end
#         # @constraint(m,Pg[gen] - Pga[gen] + unable_to_ramp_down[gen]>= -ramp*ref[:gen][gen]["pmax"]*ug[gen] - ref[:gen][gen]["pmax"]*(1 - ug[gen]))
#         # @constraint(m,Pg[gen] - Pga[gen] >= -ramp*ref[:gen][gen]["pmax"])
#         # @constraint(m,Pg[gen] - Pga[gen] <= ramp*ref[:gen][gen]["pmax"])
#     end
#     for gen in ess
#         @constraint(m,Pes[gen]<=ref[:battery][gen]["pmax_out"])
#         @constraint(m,Pes[gen]>=-ref[:battery][gen]["pmax_in"])
#         @constraint(m,Ses[gen] == Sesa[gen] - Pes[gen])
#         @constraint(m,Sesa[gen] - Pes[gen]*ref[:battery][gen]["η_in"] <= ref[:battery][gen]["SOC_max"]*ref[:battery][gen]["C_max"])
#         @constraint(m,Sesa[gen] - Pes[gen]/ref[:battery][gen]["η_out"] >= ref[:battery][gen]["SOC_min"]*ref[:battery][gen]["C_max"])
#         @constraint(m,t[gen] >= Pes[gen])
#         @constraint(m,t[gen] >= - Pes[gen])
#         @constraint(m,t[gen] <= max(ref[:battery][gen]["pmax_out"],ref[:battery][gen]["pmax_in"]))
#         @constraint(m,t[gen] >= 0)
#         @constraint(m,SocGuide[gen] >= 0)
#         @constraint(m,SocGuide[gen] <= 5*ref[:battery][gen]["C_max"]*γ_es)
#     end
#     for gen in wfs
#         @constraint(m,Pw[gen] <= Pw_max[gen])
#         @constraint(m,Pw[gen] >= 0)
#     end
#     #energy balance
#     # p_inj = Dict()
#     # for bus in keys(ref[:bus])
#     #     p = 0
#     #     if !isempty(ref[:bus_gens][bus])
#     #         p = p + sum(Pg[gen] for gen in ref[:bus_gens][bus])
#     #     end
#     #     if !isempty(ref[:bus_loads])
#     #         p = p #TODO fix bus load
#     #     end
#     # end
#     @constraint(m,data[:load] - loadCut == sum(Pw) + sum(Pg) + sum(Pes))
#     @constraint(m,loadCut>=-data[:load])
#     @constraint(m,loadCut<=data[:load])
#     @constraint(m,unbalance >= -loadCut)
#     @constraint(m,unbalance >= loadCut)
#     # @constraint(m,unbalance >= 0)
#     @constraint(m,unbalance <= data[:lo  ad])
#     return blm
# end
function makePTDF(ref)
    bus_ext = collect(keys(ref[:bus]))
    branch_ext = collect(keys(ref[:branch]))
    bus_int = 1:length(bus_ext)
    ext2int = Dict([(bus_ext[i],i) for i in bus_int])
    B = zeros(length(bus_int),length(bus_int))
    # A = zeros(length(branch_ext),length(bus_int))
    for line in branch_ext
        f_bus = ref[:branch][line]["f_bus"]
        t_bus = ref[:branch][line]["t_bus"]
        B[ext2int[f_bus],ext2int[t_bus]] = - 1/ref[:branch][line]["br_x"]
        B[ext2int[t_bus],ext2int[f_bus]] = - 1/ref[:branch][line]["br_x"]
        B[ext2int[f_bus],ext2int[f_bus]] += 1/ref[:branch][line]["br_x"]
        B[ext2int[t_bus],ext2int[t_bus]] += 1/ref[:branch][line]["br_x"]
    end
    
    # B[end,:] = sum(B[1:end-1,:],dims=1)
    Binv = pinv(B)
    PTDF = Dict()
    function cut(x)
        if abs(x)<=1e-3
            return 0.0
        else
            return x
        end
    end
    for line in branch_ext
        f_bus = ref[:branch][line]["f_bus"]
        t_bus = ref[:branch][line]["t_bus"]
        coef = map(x->cut((Binv[ext2int[f_bus],ext2int[x]] - Binv[ext2int[t_bus],ext2int[x]])/ref[:branch][line]["br_x"]),bus_ext)
        PTDF[line] = Dict(bus_ext .=> coef)
        # for bus in bus_ext
        #     PTDF[line][bus] = (Binv[ext2int[f_bus],ext2int[bus]] - Binv[ext2int[t_bus],ext2int[bus]])/ref[:branch][line]["br_x"]
        # end
    end
    ref[:PTDF] = PTDF
end
function intradayProblem(ref,data,config) #used to calculate the dual
    # m = JuMP.Model(with_optimizer(Gurobi.Optimizer))
    T,ramp,UT,DT,γ_load,γ_wind,γ_es,α,K_seg = config["T"],config["ramp"],config["UT"],config["DT"],config["γ_load"],config["γ_wind"],config["γ_es"],config["α"],config["K_seg"]
    PenaltyFactor,has_pf = config["PenaltyFactor"],config["has_pf"]
    m = JuMP.direct_model(Gurobi.Optimizer())
    aggr = []
    gens_param_lost = [gen for gen in keys(ref[:gen]) if length(ref[:gen][gen]["cost"])<2]
    gens = setdiff([gen for gen in keys(ref[:gen])],aggr)
    gens = setdiff(gens,gens_param_lost)
    ess = [gen for gen in keys(ref[:battery])]
    wfs = [gen for gen in keys(ref[:windfarm])]
    @variable(m,Pw_max[wfs],lower_bound=0)
    @variable(m,Pw_err[wfs])
    @variable(m,Pes[ess])
    @variable(m,t[ess])
    @variable(m,Ses[ess])
    @variable(m,Sesa[ess])#ancestor
    # @variable(m,SocGuide[ess])
    @variable(m,Pg[gens])
    # @variable(m,Pgu[gens])#dayahead decision
    # @variable(m,Pgl[gens])#dayahead decision
    @variable(m,Pgcost[gens],lower_bound=0)
    @variable(m,ug[gens])
    @variable(m,Pga[gens])#ancestor
    @variable(m,Pw[wfs],lower_bound=0)
    @variable(m,cost_to_go,lower_bound=0)
    @variable(m,ess_tail_cost,lower_bound=0)
    if has_pf
        @variable(m,bus_load_cut_aux[keys(ref[:load])],lower_bound=0)
        @variable(m,bus_load_cut[keys(ref[:load])])
        # @variable(m,Pf[keys(ref[:branch])])
        # @variable(m,θ[keys(ref[:bus])])
        # ref_bus = collect(ref[:ref_buses])[1][1]
        # @constraint(m,θ[ref_bus] == 0)
    else
        @variable(m,loadCut)
        @variable(m,loadCut_aux,lower_bound=0)    
    end
    #objective value
    fuel_cost = @expression(m,sum(Pgcost))
    # deg_cost = @expression(m,sum(t[gen]*ref[:battery][gen]["deg_cost"] for gen in ess))
    wind_Cur_cost = @expression(m,24/T * γ_wind*(sum(Pw_max) - sum(Pw)))
    if has_pf
        load_Cut_cost = @expression(m,24/T * γ_load*sum(bus_load_cut_aux))
    else
        load_Cut_cost = @expression(m,24/T * γ_load*loadCut_aux)
    end
    ess_cost = @expression(m,24/T *sum(t))
    obj = fuel_cost + wind_Cur_cost + load_Cut_cost + cost_to_go + ess_tail_cost + ess_cost
    m[:cost_now] = @expression(m,obj-cost_to_go)
    m[:gens] = gens
    @objective(m,Min,obj)
    for gen in wfs
        # @constraint(m,Pw_err[gen] >= -(α+0.001)*data[:wind_power][gen])
        # @constraint(m,Pw_err[gen] <= (α+0.001)*data[:wind_power][gen])
        @constraint(m,Pw_max[gen] == Pw_err[gen] + data[:wind_power][gen])
    end
    for gen in gens
        if length(ref[:gen][gen]["cost"]) == 3
            for k = 1:K_seg
                pk = ref[:gen][gen]["pmax"]*k/K_seg
                costk = (ref[:gen][gen]["cost"][2])*pk + (ref[:gen][gen]["cost"][1])*pk^2
                ratek = 2*ref[:gen][gen]["cost"][1]*pk
                @constraint(m,Pgcost[gen] >= 24/T * (costk + ratek*(Pg[gen] - pk)))
            end
        elseif length(ref[:gen][gen]["cost"]) == 2
            @constraint(m,Pgcost[gen] == 24/T * ref[:gen][gen]["cost"][1]*Pg[gen])
        else
            error("cost term must be 2 or 3")
        end
        # @constraint(m,Pg[gen] <= Pgu[gen])
        # @constraint(m,Pg[gen] >= Pgl[gen])
        @constraint(m,Pg[gen] <= ug[gen] * ref[:gen][gen]["pmax"])
        @constraint(m,Pg[gen] >= 0)
        @constraint(m,Pg[gen] - Pga[gen] >= -24/T * ramp*ref[:gen][gen]["pmax"]*ug[gen] - ref[:gen][gen]["pmax"]*(1 - ug[gen]))
        @constraint(m,Pg[gen] - Pga[gen] <= 24/T * ramp*ref[:gen][gen]["pmax"])
    end
    for gen in ess
        @constraint(m,Pes[gen]<=ref[:battery][gen]["pmax_out"])
        @constraint(m,Pes[gen]>=-ref[:battery][gen]["pmax_in"])
        @constraint(m,Ses[gen] == Sesa[gen] - 24/T * Pes[gen])
        @constraint(m,Sesa[gen] - 24/T * Pes[gen]*ref[:battery][gen]["η_in"] <= ref[:battery][gen]["SOC_max"]*ref[:battery][gen]["C_max"])
        @constraint(m,Sesa[gen] - 24/T * Pes[gen]/ref[:battery][gen]["η_out"] >= ref[:battery][gen]["SOC_min"]*ref[:battery][gen]["C_max"])
        @constraint(m,t[gen]>= ref[:battery][gen]["deg_cost"]*Pes[gen])
        @constraint(m,t[gen]>= - ref[:battery][gen]["deg_cost"]*Pes[gen])
        # @constraint(m,SocGuide[gen] >= 0)
    end
    for gen in wfs
        @constraint(m,Pw[gen] <= Pw_max[gen])
    end
    #energy balance
    if has_pf
                #     transmission constraint 
        #     #energy balance
        total_load = length(ref[:load])
        p_inj = Dict()
        for bus in keys(ref[:bus])
            p = AffExpr()
            if !isempty(ref[:bus_gens][bus])
                add_to_expression!(p,sum(Pg[gen] for gen in ref[:bus_gens][bus])) 
            end
            p_inj[bus] = p
        end
        for load in keys(ref[:load])
            info = ref[:load][load]
            add_to_expression!(p_inj[info["load_bus"]], - data[:load]/total_load + bus_load_cut[load])
            @constraint(m,bus_load_cut[load]<=data[:load]/total_load)
            @constraint(m,bus_load_cut[load]>= - data[:load]/total_load)
            @constraint(m,bus_load_cut_aux[load] >= bus_load_cut[load])
            @constraint(m,bus_load_cut_aux[load] >= - bus_load_cut[load])
        end
        for gen in wfs
            idx = ref[:windfarm][gen]
            add_to_expression!(p_inj[idx], Pw[gen])
        end
        for gen in ess
            idx = ref[:battery_location][gen]
            add_to_expression!(p_inj[idx], Pes[gen])
        end
        for line in keys(ref[:branch])
            f_bus = ref[:branch][line]["f_bus"]
            t_bus = ref[:branch][line]["t_bus"]
            br_x = ref[:branch][line]["br_x"]
            # @constraint(m,Pf[line] == (θ[f_bus]-θ[t_bus])/br_x)
            if haskey(ref[:branch][line],"rate_a") # Which means power flow limits are enforced
                pf = sum(v for (k,v) in ref[:PTDF][line])
                # pf1 = reduce(+,pf)
                @constraint(m,pf<=ref[:branch][line]["rate_a"])
                @constraint(m,pf>=-ref[:branch][line]["rate_a"])
            end
        end
        @constraint(m,sum(p_inj[k] for k in keys(p_inj)) == 0)
        # for bus in keys(ref[:bus])
        #     @constraint(m,p_inj[bus] == 0)
        # end
        m[:P_inj] = p_inj
    else
        @constraint(m,data[:load] - loadCut == sum(Pw) + sum(Pg) + sum(Pes))
        @constraint(m,loadCut <= data[:load])
        @constraint(m,loadCut >= -data[:load])
        @constraint(m,loadCut_aux >= loadCut)
        @constraint(m,loadCut_aux >= -loadCut)    
    end
    # @constraint(m,upper_bound,cost_to_go >= 0)
    # construct the overestimator
    mo = JuMP.Model(with_optimizer(Gurobi.Optimizer))
    @variable(mo,μ)
    @variable(mo,τ_Pg[gens],upper_bound=2*γ_load,lower_bound=-2*γ_load)
    @variable(mo,τ_Ses[ess],upper_bound=2*γ_load,lower_bound=-2*γ_load)
    m[:states] = [Pg[:].data;Ses[:].data]#TODO extend this method
    m[:states_copy] = [Pga[:].data;Sesa[:].data]#TODO extend this method
    # m[:sum_yk] = - 1
    @variable(m,τu[1:length(m[:states])],lower_bound=0)
    @variable(m,τl[1:length(m[:states])],lower_bound=0)
    @constraint(m,sum_states[i=1:length(m[:states])],m[:states][i] + τu[i] - τl[i] == 0)
    # m[:sum_states] = [m[:states][i] + τu[i] - τl[i] for i in 1:length(m[:states])]
    # m[:sum_states_cons] = []
    @constraint(m,upper_bound,cost_to_go >= sum(τu[i]*γ_load*PenaltyFactor for i in 1:length(m[:states])) + sum(τl[i]*γ_load*PenaltyFactor for i in 1:length(m[:states])))
    # @objective(mo,Max,μ+sum(τ_Pg)+sum(τ_Ses))
    m[:cost_to_go_now] = cost_to_go
    m[:overestimator] = mo
    m[:converging] = false
    m[:penalty] = [1000.0 for i in 1:length(m[:states])]
    m[:α] = config["α"]
    return m
end
function PolygonUncertaintySet(center,covariance,Γ,s=3)
    nx = length(center)
    @assert 1<=Γ<=nx
    if nx == 1
        vertice = [-s*covariance,s*covariance]
        return vertice
    else
        m = Model()
        @variable(m,-1<=x[1:nx]<=1)
        for idx in 1:2^nx
            indicator = bitstring(idx)[end-nx+1:end]
            indicator = [parse(Int,indicator[i]) * 2 - 1 for i in 1:nx]
            @constraint(m,sum(x[i]*indicator[i] for i in 1:nx)<=Γ)
        end
        poly = polyhedron(m,CDDLib.Library())
        vertice = []
        for v in points(poly)
            push!(vertice,center+s*covariance^0.5*v)
            # push!(vertice,v)
        end
        return vertice
    end
end
function getUnionVertice(unc::Array{})
    vall = []
    for vertice in unc
        vall = union(vall,vertice)
    end
    if length(vall[1]) == 1
        return [[x] for x in vall]
    else
        P = vrep([v for v in vall])
        Q = polyhedron(P)
        removevredundancy!(Q)
        return Q.vrep.points.points
    end
end
function dayaheadProblem(ref,data,config)
    @info("constructing dayahead model")
    T,ramp,UT,DT,γ_load,γ_wind,γ_es,α,K_seg = config["T"],config["ramp"],config["UT"],config["DT"],config["γ_load"],config["γ_wind"],config["γ_es"],config["α"],config["K_seg"]
    PenaltyFactor = config["PenaltyFactor"]
    m = JuMP.Model(with_optimizer(Gurobi.Optimizer,MIPGap = 0.01))
    aggr = []
    gens_param_lost = [gen for gen in keys(ref[:gen]) if length(ref[:gen][gen]["cost"])<2]
    gens = setdiff([gen for gen in keys(ref[:gen])],aggr)
    gens = setdiff(gens,gens_param_lost)
    pmax_dict = [(gen,ref[:gen][gen]["pmax"]) for gen in gens]
    gen_top_sort_by_pmax = Base.sort(pmax_dict,by=x->x[2],rev=true)[1:min(convert(Int,floor((length(gens)/2))),30)]
    gen_top_sort_by_pmax = [x[1] for x in gen_top_sort_by_pmax]
    m[:gens] = gens
    ess = [gen for gen in keys(ref[:battery])]
    wfs = [gen for gen in keys(ref[:windfarm])]
    T = length(data)
    Ts = range(1,stop=length(data))
    @variable(m,Pes[ess,Ts],lower_bound=0)
    @variable(m,t[ess,Ts])
    @variable(m,Ses[ess,Ts])
    @variable(m,Pg[gens,Ts])
    @variable(m,Pgcost[gens,Ts],lower_bound=0)
    @variable(m,Pgl[gens,Ts],lower_bound=0)
    @variable(m,Pgu[gens,Ts],lower_bound=0)
    @variable(m,ug[gens,Ts],Bin)
    @variable(m,Pw[wfs,Ts],lower_bound=0)
    @variable(m,loadCut[Ts],lower_bound=0)
    @variable(m,cost_to_go,lower_bound=0)
    @variable(m,ug_aux[gens,Ts])
    @constraint(m,PguLimit[gen in gens,t in 1:T],Pgu[gen,t] <= ug[gen,t]*ref[:gen][gen]["pmax"])
    @constraint(m,PglLimit[gen in gens,t in 1:T],Pgl[gen,t] >= ug[gen,t]*ref[:gen][gen]["pmin"])
    for t in Ts
        for gen in gens
            for k = 1:K_seg
                pk = ref[:gen][gen]["pmax"]*k/K_seg
                costk = (ref[:gen][gen]["cost"][2])*pk + (ref[:gen][gen]["cost"][1])*pk^2
                ratek = 2*ref[:gen][gen]["cost"][1]*pk
                @constraint(m,Pgcost[gen,t] >= costk + ratek*(Pg[gen,t] - pk))
            end
            @constraint(m,Pg[gen,t] >= Pgl[gen,t])
            @constraint(m,Pg[gen,t] <= Pgu[gen,t])
            if t >= 2
                #RAMP COSTRAINT
                @constraint(m,Pgl[gen,t] - Pgu[gen,t-1] >= -24/T * ramp*ref[:gen][gen]["pmax"]*ug[gen] - ref[:gen][gen]["pmax"]*(1 - ug[gen]))
                @constraint(m,Pgu[gen,t] - Pgl[gen,t-1] <= 24/T * ramp*ref[:gen][gen]["pmax"])
                #STARTUP-SHUTDOWN CONSTRAINT
                if t <= T - UT + 1
                    @constraint(m,sum(ug[gen,k] for k in t:t+UT-1) >= T/24 * UT*(ug[gen,t] - ug[gen,t-1]))
                else
                    @constraint(m,sum(ug[gen,k] - (ug[gen,t] - ug[gen,t-1]) for k in t:T) >= 0)
                end
                if t <= T - DT + 1
                    @constraint(m,sum(1 - ug[gen,k] for k in t:t+DT-1) >= T/24 * DT*(ug[gen,t-1] - ug[gen,t]))
                else
                    @constraint(m,sum(1 - ug[gen,k] - (ug[gen,t-1] - ug[gen,t]) for k in t:T) >= 0)
                end
                @constraint(m,ug_aux[gen,t] >= ug[gen,t] - ug[gen,t-1])
                @constraint(m,ug_aux[gen,t] >= 0)
            end
        end
        for gen in ess
            @constraint(m,Pes[gen,t]<=ref[:battery][gen]["pmax_out"])
            @constraint(m,Pes[gen,t]>=-ref[:battery][gen]["pmax_in"])
            if t == 1
                @constraint(m,Ses[gen,t] == ref[:battery][gen]["SOC_int"]*ref[:battery][gen]["C_max"])
            else
                @constraint(m,Ses[gen,t] == (1 - ref[:battery][gen]["sigma"] * 24/T)*Ses[gen,t-1] - 24/T * Pes[gen,t-1])
            end
            if t == T
                @constraint(m,Ses[gen,t] == Ses[gen,1])
            end
            @constraint(m,Ses[gen,t] <= ref[:battery][gen]["SOC_max"]*ref[:battery][gen]["C_max"])
            @constraint(m,Ses[gen,t] >= ref[:battery][gen]["SOC_min"]*ref[:battery][gen]["C_max"])
        end
        for gen in wfs
            @constraint(m,Pw[gen,t] <= data[t][:wind_power][gen])
        end
        #power balance
        @constraint(m,loadCut[t] <= data[t][:load])
        @constraint(m,data[t][:load] - loadCut[t]== sum(Pw[:,t]) + sum(Pg[:,t]) + sum(Pes[:,t]))
    end
    startup_cost = @expression(m,sum(ref[:gen][gen]["startup"]*ug_aux[gen,t] for gen in gens,t in 2:T))
    fuel_cost = @expression(m,sum(Pgcost))
    run_cost = 0
    interval_reward = @expression(m,sum(Pgu[gen,t] - Pgl[gen,t] for gen in gens, t in Ts))
    # try
    #     run_cost = @expression(m,sum((ref[:gen][gen]["cost"][3]) * ug[gen,t] for gen in gens,t in Ts))# May slow down the program
    # catch exception
    #     @info("runtime cost not available, set to zero")
    # end
    wind_total = sum(data[t][:wind_power][gen] for t in Ts,gen in wfs)
    wind_Cur_cost = @expression(m,γ_wind*(wind_total - sum(Pw)))
    m[:wind_Cur_cost] = wind_Cur_cost
    m[:cost_now] = startup_cost + run_cost
    load_Cut_cost = @expression(m,γ_load*sum(loadCut))
    @constraint(m,cost_to_go >= fuel_cost + wind_Cur_cost + load_Cut_cost)
    @objective(m,Min,startup_cost + run_cost - interval_reward + cost_to_go)
    return m
end

function prior_list_modification(dayahead,ref)
    T = length(dayahead[:loadCut])
    # @suppress_out begin
    model = deepcopy(dayahead)
    for gen in model[:gens]
        for t in 1:T
            fix(model[:ug][gen,t],1)
        end
    end
    for gen in model[:gens]
        for t in 1:T
            set_normalized_coefficient(model[:PglLimit][gen,t],model[:ug][gen,t],0)
        end
    end
    optimize!(model)
    num_have_to_be_on = 0
    num_have_to_be_off = 0
    for gen in model[:gens]
        on_flag = true
        off_flag = true
        for t in 1:T
            if value(model[:Pg][gen,t]) < 0.65 * ref[:gen][gen]["pmax"]
                on_flag = false
                break
            end
        end
        for t in 1:T
            if value(model[:Pg][gen,t]) > 0.1 * ref[:gen][gen]["pmax"]
                off_flag = false
                break
            end
        end
        if on_flag
            for t in 1:T
                fix(dayahead[:ug][gen,t],1)
            end
            num_have_to_be_on += 1
            # @info("fixed on")
        elseif off_flag
            num_have_to_be_off += 1
            for t in 1:T
                fix(dayahead[:ug][gen,t],0)
            end
            # @info("fixed off")
        else
            nothing
        end
    end
    @info("The number of always-on units is $num_have_to_be_on")
    @info("The number of always-off units is $num_have_to_be_off")
    # end
    # return dayahead
end

function fix_tail()
    model = Main.RTD.intraday[end]
    tol = 0.2
    ref = Main.RTD.case_dict
    ess = keys(ref[:battery])
    @variable(model,Qes[ess])
    @constraint(model,Qes_lower[gen = ess],Qes[gen] >= -ref[:battery][gen]["penalty"] * (model[:Ses][gen] -  (1-tol)*ref[:battery][gen]["SOC_int"] * ref[:battery][gen]["C_max"]))
    @constraint(model,Qes_upper[gen = ess],Qes[gen] >= ref[:battery][gen]["penalty"] * (model[:Ses][gen] - (1+tol)*ref[:battery][gen]["SOC_int"] * ref[:battery][gen]["C_max"]))
    @constraint(model,model[:ess_tail_cost] >= sum(Qes))
    model = Main.RTD.intradayMax[end]
    # @constraint(model,Qes_equal[gen = ess],model[:Ses][gen] ==  ref[:battery][gen]["SOC_int"] * ref[:battery][gen]["C_max"])
    @variable(model,Qes[ess])
    @constraint(model,Qes_lower[gen = ess],Qes[gen] >= -ref[:battery][gen]["penalty"] * (model[:Ses][gen] -  (1-tol)*ref[:battery][gen]["SOC_int"] * ref[:battery][gen]["C_max"]))
    @constraint(model,Qes_upper[gen = ess],Qes[gen] >= ref[:battery][gen]["penalty"] * (model[:Ses][gen] - (1+tol)*ref[:battery][gen]["SOC_int"] * ref[:battery][gen]["C_max"]))
    @constraint(model,model[:ess_tail_cost] >= sum(Qes))
    return model
end
function fixall(t)
    fixall(Main.RTD,t)
end
function fixall(model::RealTimeDispatchModel,t::Int)
    T = length(model.intraday)
    for gen in model.dayahead[:gens]
        # fix(case_base[:Pgu][gen], value(dayahead[:Pgu][gen,t])) # box method
        # fix(case_base[:Pgl][gen], value(dayahead[:Pgl][gen,t]))
        # fix(case_max[:Pgu][gen], value(dayahead[:Pgu][gen,t]))
        # fix(case_max[:Pgl][gen], value(dayahead[:Pgl][gen,t]))
        # x = model.dayahead[:ug][gen,t]
        x = 1
        fix(model.intraday[t][:ug][gen],x)
        fix(model.intradayMax[t][:ug][gen],x)
    end
    # fix the current state varible
    for gen in model.dayahead[:gens]
        if t != 1
            x = value(model.intraday[t-1][:Pg][gen])
            fix(model.intraday[t][:Pga][gen],x)
            fix(model.intradayMax[t][:Pga][gen],x)
        end
    end
    for gen in keys(model.case_dict[:battery])
        if t != 1
            Sesa_real =  (1 - model.case_dict[:battery][gen]["sigma"] * 24/T) * value(model.intraday[t-1][:Sesa][gen]) -
            24/T * min(value(model.intraday[t-1][:Pes][gen])*model.case_dict[:battery][gen]["η_in"],0) -
            24/T * max(value(model.intraday[t-1][:Pes][gen])/model.case_dict[:battery][gen]["η_out"],0)
            # @constraint(case_min,case_min[:Sesa][gen] == Sesa_real)
            fix(model.intraday[t][:Sesa][gen],Sesa_real)
            fix(model.intradayMax[t][:Sesa][gen],Sesa_real)
            # @everywhere fix(Main.intradayMax[$t][:Sesa][$gen],$Sesa_real)
        else
            soc_int = model.case_dict[:battery][gen]["SOC_int"]*model.case_dict[:battery][gen]["C_max"]
            fix(model.intraday[t][:Sesa][gen],soc_int)
            fix(model.intradayMax[t][:Sesa][gen],soc_int)
            # @everywhere fix(Main.intradayMax[$t][:Sesa][$gen],$soc_int)
        end
    end
end
function eval_wst_case(t,wst_vertex)
    eval_wst_case(Main.RTD,t,wst_vertex)
    return 0
end
function eval_wst_case(model::RealTimeDispatchModel,t,wst_vertex)
    case_base = model.intraday[t]
    case_max = model.intradayMax[t]
    for (idx,gen) in enumerate(keys(model.case_dict[:windfarm]))
        # @info(wst_vertex)
        fix(case_base[:Pw_err][gen],wst_vertex[idx]*case_base[:α]*model.data[t][:wind_power][gen])
        fix(case_max[:Pw_err][gen],wst_vertex[idx]*case_base[:α]*model.data[t][:wind_power][gen])
    end
    optimize!(case_base)
    optimize!(case_max)
    @assert termination_status(case_base) == MOI.OPTIMAL
    @assert termination_status(case_max) == MOI.OPTIMAL 
    return 0
end
function fix_and_optimize(t,i)
    case_tmp = Main.RTD.intradayMax[t] # TODO: not sure
    # set_opt【imizer(case_tmp,with_optimizer(Gurobi.Optimizer,Main.env))
    vertex = Main.RTD.vertice[t][i]
    windpower = Main.RTD.data[t][:wind_power]
    for (idx,gen) in enumerate(keys(Main.RTD.case_dict[:windfarm]))
        @assert abs(vertex[idx])*Problems.α <= 0.5
        fix(case_tmp[:Pw_err][gen],vertex[idx]*Problems.α*windpower[gen])
    end
    @suppress optimize!(case_tmp)
    # @assert termination_status(case_tmp) == MOI.OPTIMAL
    if objective_value(case_tmp) <= 0
        @warn("objective less than zero with vertex $vertex")
    end
    # if objective_value(case_max) >= wst_case_value
    #     wst_case_value = objective_value(case_max)
    #     wst_vertex = vertex
    #     # wst_case_max = deepcopy(case_max)
    # end
    return objective_value(case_tmp)
    # @info("$i finished in $(round(t2-t1; digits=2)) seconds on worker $(myid())")

end
function fix_and_optimize(model::RealTimeDispatchModel,t,i)
    case_tmp = model.intradayMax[t] # TODO: not sure
    # set_opt【imizer(case_tmp,with_optimizer(Gurobi.Optimizer,Main.env))
    vertex = model.vertice[t][i]
    windpower = model.data[t][:wind_power]
    for (idx,gen) in enumerate(keys(model.case_dict[:windfarm]))
        @assert abs(vertex[idx])*case_tmp[:α] <= 0.5
        fix(case_tmp[:Pw_err][gen],vertex[idx]*case_tmp[:α]*windpower[gen])
        # try
        #     fix(case_tmp[:Pw_err][gen],vertex[idx]*α*data[:wind_power][gen])
        # catch e
        #     gurobi_env = Gurobi.Env()
        #     set_optimizer(case_tmp,with_optimizer(Gurobi.Optimizer,gurobi_env))
        #     fix(case_tmp[:Pw_err][gen],vertex[idx]*α*data[:wind_power][gen])
        # end
    end
    @suppress optimize!(case_tmp)
    # @assert termination_status(case_tmp) == MOI.OPTIMAL
    if objective_value(case_tmp) <= 0
        @warn("objective less than zero with vertex $vertex")
    end
    # if objective_value(case_max) >= wst_case_value
    #     wst_case_value = objective_value(case_max)
    #     wst_vertex = vertex
    #     # wst_case_max = deepcopy(case_max)
    # end
    return objective_value(case_tmp)
    # @info("$i finished in $(round(t2-t1; digits=2)) seconds on worker $(myid())")

end
function add_upper_bound(t,n_iter)
    add_upper_bound(Main.RTD,t,n_iter)
    return 0
end
function add_upper_bound(model::RealTimeDispatchModel,t,n_iter)
    dual_of_obj = @variable(model.intradayMax[t],lower_bound=0)
    if n_iter == 1
        @constraint(model.intradayMax[t],sum_y,dual_of_obj == 1)
    end
    v_upper = objective_value(model.intradayMax[t+1])
    #modify objective
    set_normalized_coefficient(model.intradayMax[t][:upper_bound],dual_of_obj,-v_upper)
    # modify constraint  ∑y_k == 1
    set_normalized_coefficient(model.intradayMax[t][:sum_y],dual_of_obj,1)
    # add constraint
    # x_i - ∑y_k*x_{ki} + τu_i - τl_i == 0 for x_i in states variables
    for i in 1:length(model.intradayMax[t][:sum_states])
        set_normalized_coefficient(model.intradayMax[t][:sum_states][i],dual_of_obj, - value(model.intraday[t][:states][i]))
    end
end
function add_lower_bound(t)
    add_lower_bound(Main.RTD,t)
    return 0
end
function add_lower_bound(model::RealTimeDispatchModel,t)
    lower_cut = AffExpr()
    for i in 1:length(model.intraday[t][:states])
        pi = dual(FixRef(model.intraday[t+1][:states_copy][i]))
        model.intraday[t][:penalty][i] = max(model.intraday[t][:penalty][i],1.01*abs(pi))
        add_to_expression!(lower_cut,pi * (model.intraday[t][:states][i] - value(model.intraday[t][:states][i])))
    end
    # for gen in keys(model.case_dict[:battery])
    #     #firstly the optimistic(base) case
    #     pi =  dual(FixRef(model.intraday[t+1][:Sesa][gen]))#TODO 
    #     model.intraday[t][:penalty] = max(model.intraday[t][:penalty],1.01*abs(pi))
    #     add_to_expression!(lower_cut,pi * (model.intraday[t][:Ses][gen] - value(model.intraday[t][:Ses][gen])))
    # end
    # for gen in model.intraday[t][:gens]
    #     pi = dual(FixRef(model.intraday[t+1][:Pga][gen]))
    #     model.intraday[t][:penalty] = max(model.intraday[t][:penalty],1.01*abs(pi))
    #     add_to_expression!(lower_cut,pi * (model.intraday[t][:Pg][gen] - value(model.intraday[t][:Pg][gen])))
    # end
    # v_lower = sum(value(model.intraday[τ][:cost_now]) for τ in t+1:T)
    v_lower = objective_value(model.intraday[t+1])
    if v_lower <= objective_value(model.intradayMax[t+1])
        global number_of_reject
        number_of_reject = 0 
        add_to_expression!(lower_cut,v_lower)
        @constraint(model.intraday[t],model.intraday[t][:cost_to_go] >= lower_cut)
    else
        # global number_of_reject
        # number_of_reject += 1
        if model.intraday[t][:converging]
            # add_to_expression!(lower_cut,objective_value(model.intradayMax[t+1]))
            add_to_expression!(lower_cut,v_lower)
            @constraint(model.intraday[t],model.intraday[t][:cost_to_go] >= lower_cut)
        else
            model.intraday[t][:converging] = true
            # add_to_expression!(lower_cut,objective_value(model.intradayMax[t+1]))
            add_to_expression!(lower_cut,v_lower)
            @constraint(model.intraday[t],model.intraday[t][:cost_to_go] >= lower_cut)
        end
    end
end
function vertex_enumeration_primal(model::RealTimeDispatchModel,t)
    start = time()
    wst_case_value = 0
    wst_vertex = 0
    wst_case_max = Nothing
    vertice = model.vertice[t]
    if nprocs == 1
        objectives_along_vertice = zeros(length(vertice))
    else
        objectives_along_vertice = SharedArray{Float64}(length(vertice))
    end
    # m = deepcopy(case_max)
    # model_shared = [m for t in 1:length(vertice)]
    if nprocs() == 1
        for i in 1:length(vertice) #vertex enumeration getting the worst case
            val = fix_and_optimize(model,t,i)
            objectives_along_vertice[i] = val
        end
    else
        @sync @distributed for i in 1:length(vertice) #vertex enumeration getting the worst case
            val = Main.Problems.fix_and_optimize(t,i)
            objectives_along_vertice[i] = val
        end
    end
    ending = time()
    #@info(ending - start) 
    global total_solves 
    total_solves += length(vertice)
    ind = findmax(objectives_along_vertice)[2]
    wst_vertex = vertice[ind]
    # if t == 8
    #     top_vertex_index = sort([x for x in 1:length(vertice)],by=x->objectives_along_vertice[x],rev=true)
    #     @info(top_vertex_index)
    # end
    return wst_vertex
end
function vertex_enumeration_dual(case_max,vertice,case_dict,data,t,N_ITER)
    wst_case_value = 0
    wst_vertex = 0
    for vertex in vertice #vertex enumeration getting the worst case
        for (idx,gen) in enumerate(keys(case_dict[:windfarm]))
            @assert abs(vertex[idx])*α <= 0.5
            fix(case_max[:Pw_err][gen],vertex[idx]*α*data[:wind_power][gen])
        end
        # benders decomposition for solving the worst scenario Problems,cuts can be shared among vertice
        UB = 1e9
        LB = 1
        n_benders = 0
        while true
            n_benders += 1
            @assert n_benders <= 500
            optimize!(case_max)#get the worst scenario
            global total_solves
            total_solves += 1
            if termination_status(case_max) != MOI.OPTIMAL
                error(termination_status(case_max))
            end
            if t == T
                break
            end
            if abs(UB - LB)/LB <= 1e-3
                # println("Benders iteration number = $n_benders")
                # println("GAP = $(UB-LB)")
                break
            end
            LB = max(LB,objective_value(case_max))
            mo = case_max[:overestimator]
            # for gen in keys(mo[:τ_Pg])
            #     set_objective_coefficient(mo,mo[:τ_Pg][gen],value(case_max[:Pg][gen]))
            # end
            # for gen in keys(mo[:τ_Ses])
            #     set_objective_coefficient(mo,mo[:τ_Ses][gen],value(case_max[:Ses][gen]))
            # end
            @objective(mo,Max,mo[:μ] +
            sum(mo[:τ_Pg][gen]*value(case_max[:Pg][gen]) for gen in keys(case_max[:Pg])) +
            sum(mo[:τ_Ses][gen]*value(case_max[:Ses][gen]) for gen in keys(case_max[:Ses])))
            optimize!(mo)
            global total_solves
            total_solves += 1
            @assert termination_status(mo) != MOI.INFEASIBLE
            if N_ITER == 1
                break
            else
                @assert termination_status(mo) == MOI.OPTIMAL
                UB = min(UB,objective_value(case_max) - value(case_max[:cost_to_go]) + objective_value(mo))
                #add cut
                @constraint(case_max,case_max[:cost_to_go] >= objective_value(mo) +
                sum((case_max[:Pg][gen] - value(case_max[:Pg][gen]))*value(mo[:τ_Pg][gen]) for gen in case_max[:gens]) +
                sum((case_max[:Ses][gen] - value(case_max[:Ses][gen]))*value(mo[:τ_Ses][gen]) for gen in keys(case_dict[:battery])))
            end
        end
        if objective_value(case_max) <= 0
            @warn objective_value(case_max)
        end
        if objective_value(case_max) >= wst_case_value
            wst_case_value = objective_value(case_max)
            wst_vertex = vertex
            # wst_case_max = deepcopy(case_max)
        end
    end
    for (idx,gen) in enumerate(keys(case_dict[:windfarm]))
        # @info(wst_vertex)
        fix(case_max[:Pw_err][gen],wst_vertex[idx]*α*data[:wind_power][gen])
    end
    optimize!(case_max)
    global total_solves
    total_solves += 1
    return case_max,wst_vertex
end
function ForwardPassPrimal(model::RealTimeDispatchModel,N_ITER,start,stop)
    additional = Dict()
    additional[:gaphourly] = []
    # additional[:upper] = []
    gap = 0
    differ = 0
    @assert length(model.intraday) == length(model.vertice)
    # intradayMax = deepcopy(intradayMax)
    # intraday = deepcopy(intraday)
    t11 = 0 #OK
    t22 = 0 #reasonable
    t33 = 0
    sceanio_not_change = true
    for t in start:stop #前推步骤
        # case_base = deepcopy(intraday[t])
        # case_max = deepcopy(intradayMax[t])
        # fix the dayahead decision
        t1 = time()
        if nprocs() == 1
            fixall(model,t)
        else
            @everywhere Problems.fixall($t) # RealTimeDispatch model are exposed to all workers in Main Module
        end
        t2 = time()
        t11 += t2 - t1
        # vertex_enumeration //makes case_max dirty
        # case_max_clean = deepcopy(case_max)
        t1 = time()
        if sceanio_not_change_num >= 10
            wst_vertex = model.intradayMax[t][:wst_vertex]
        else
            wst_vertex = vertex_enumeration_primal(model,t)
        end
        if N_ITER >= 2
            if wst_vertex != model.intradayMax[t][:wst_vertex]
                sceanio_not_change = false
            end
        end
        model.intradayMax[t][:wst_vertex] = wst_vertex
        t2 = time()
        t22 += t2 - t1
        # fix the worst vertex and calculate the intraday response
        t1 = time()
        calls = []
        if nprocs() == 1
            @suppress_out eval_wst_case(model,t,wst_vertex)
        else
            for w in procs()
                r = @spawnat w @suppress_out eval_wst_case(t,wst_vertex)
                push!(calls,r)
            end
            # @everywhere @suppress_out Problems.eval_wst_case($t,$wst_vertex)
            for (i,w) in enumerate(procs())
                fetch(calls[i])
            end
        end
        t2 = time()
        t33 += t2 - t1
        # @suppress optimize!(case_base)
        # @suppress optimize!(case_max)
        # broadcast the state vars
        global total_solves
        total_solves += 2
        gapt  = value((model.intradayMax[t][:cost_to_go]) - value(model.intraday[t][:cost_to_go]))/value(model.intradayMax[t][:cost_to_go])
        push!(additional[:gaphourly],round(gapt;digits=3))
    end
    additional[:UpperBound] = objective_value(model.intradayMax[1])
    additional[:LowerBound] = objective_value(model.intraday[1])
    # additional[:UpperBound] = value(model.intradayMax[1][:cost_to_go])
    # additional[:LowerBound] = value(model.intraday[1][:cost_to_go])
    additional[:Gap] = (additional[:UpperBound] - additional[:LowerBound])/additional[:UpperBound]
    # @info("sceanio_not_change = $sceanio_not_change")
    # @info(additional[:gaphourly])
    global sceanio_not_change_num
    if sceanio_not_change
        sceanio_not_change_num += 1
    else
        sceanio_not_change_num = 0
    end
    # @warn(t11,t22,t33)
    return additional
end
function ForwardPassDual(dayahead::JuMP.Model,intraday::Array{JuMP.Model},intradayMax::Array{JuMP.Model},vertice::Array{},case_dict,data,N_ITER)
    T = length(intraday)
    additional = Dict()
    additional[:upper] = []
    gap = 0
    differ = 0
    @assert length(intraday) == length(vertice)
    intradayMaxToken = deepcopy(intradayMax)
    # intraday = deepcopy(intraday)
    for t in 1:length(intraday) #前推步骤
        case_base = deepcopy(intraday[t])
        case_max = deepcopy(intradayMax[t])
        # fix the dayahead decision
        if t == 1
            case_base = fixall(case_base,Nothing,dayahead,case_dict,t)
            case_max = fixall(case_max,Nothing,dayahead,case_dict,t)
        else
            case_base = fixall(case_base,intraday[t-1],dayahead,case_dict,t)
            case_max = fixall(case_max,intraday[t-1],dayahead,case_dict,t)
        end
        # vertex_enumeration //makes case_max dirty
        case_max_clean = deepcopy(case_max)
        case_max,wst_vertex = vertex_enumeration_dual(case_max,vertice[t],case_dict,data[t],t,N_ITER)
        # fix the worst vertex and calculate the intraday response
        for (idx,gen) in enumerate(keys(case_dict[:windfarm]))
            # @info(wst_vertex)
            fix(case_base[:Pw_err][gen],wst_vertex[idx]*α*data[t][:wind_power][gen])
        end
        optimize!(case_base)
        global total_solves
        total_solves += 1
        @assert termination_status(case_base) == MOI.OPTIMAL
        if t == 1
            gap = (objective_value(case_max) - objective_value(case_base))/objective_value(case_max)
        end
        differ += value(case_max[:cost_now]) - value(case_base[:cost_now])
        case_max[:cost_to_go_copy] =  value(case_max[:cost_to_go])
        push!(additional[:upper],objective_value(case_max))
        intradayMax[t] = case_max_clean
        intraday[t] = case_base
        intradayMaxToken[t] = case_max
        # @info("t = $t \n    objective: $(objective_value(case_base))\n    unbalance:$(value.(case_base[:loadCut_aux]))")
        # obj_now[t] = objective_value(case_base) - value(case_base[:cost_to_go])
    end
    additional[:gap] = gap
    additional[:differ] = differ/sum(value(intraday[t][:cost_now]) for t in 1:T)
    return intraday,intradayMax,intradayMaxToken,additional
end

function BackwardPassPrimal(model::RealTimeDispatchModel,N_ITER,start=1,stop=96)
    # intraday = deepcopy(intraday)
    # intradayMax = deepcopy(intradayMax)
    for t in [stop-x+start for x in start+1:stop]#回代步骤
        # @info(t)
        # **update the overestimator**
        # case_max_clean = deepcopy(case_max)
        if sceanio_not_change_num >= 10
            wst_vertex = model.intradayMax[t+1][:wst_vertex]
        else
            wst_vertex = vertex_enumeration_primal(model,t+1)
        end
        # wst_vertex = model.intradayMax[t+1][:wst_vertex]
        if nprocs() == 1
            eval_wst_case(model,t+1,wst_vertex)
        else
            calls = []
            for w in procs()
                r = @spawnat w @suppress_out eval_wst_case(t+1,wst_vertex)
                push!(calls,r)
            end
            for (i,w) in enumerate(procs())
                fetch(calls[i])
            end
        end
        global total_solves
        total_solves += 1
        if objective_value(model.intradayMax[t+1]) < value(model.intradayMax[t][:cost_to_go]) ||  N_ITER == 1
            if nprocs() == 1
                add_upper_bound(model,t,N_ITER)
            else
                calls = [] 
                for w in procs()
                    r = @spawnat w @suppress_out add_upper_bound(t,N_ITER)
                    push!(calls,r)
                end
                for (i,w) in enumerate(procs())
                    fetch(calls[i])
                end
            end
        end
        # **solve the updated lower problem**
        # for (idx,gen) in enumerate(keys(model.case_dict[:windfarm]))
        #     fix(model.intraday[t+1][:Pw_err][gen],wst_vertex[idx]*α*model.data[t+1][:wind_power][gen])
        # end
        # optimize!(model.intraday[t+1])
        # global total_solves
        # total_solves += 1
        # **update the underestimator**
        if objective_value(model.intraday[t+1]) > value(model.intraday[t][:cost_to_go]) ||  N_ITER == 1
            if nprocs() == 1
                add_lower_bound(model,t)
            else
                calls = []
                for w in procs()
                    r = @spawnat w @suppress_out add_lower_bound(t)
                    push!(calls,r)
                end
                for (i,w) in enumerate(procs())
                    fetch(calls[i])
                end
            end
        end
    end
    for t in start:stop
        for i in 1:length(model.intradayMax[t][:states]) 
            set_normalized_coefficient(model.intradayMax[t][:upper_bound],model.intradayMax[t][:τu][i],-model.intraday[t][:penalty][i])
            set_normalized_coefficient(model.intradayMax[t][:upper_bound],model.intradayMax[t][:τl][i],-model.intraday[t][:penalty][i])
        end
    end
end

function BackwardPassDual(intraday::Array{JuMP.Model},intradayMax::Array{JuMP.Model},intradayMaxToken::Array{JuMP.Model},vertice::Array{},case_dict,data,N_ITER,start=1,stop=T)
    # intraday = deepcopy(intraday)
    # intradayMax = deepcopy(intradayMax)
    T = length(intraday)
    for t in [T-x+1 for x in 2:T]#回代步骤
        # @info(t)
        # **update the overestimator**
        case_max = intradayMax[t+1]
        case_max_clean = deepcopy(case_max)
        case_max,wst_vertex = vertex_enumeration_dual(case_max,vertice[t+1],case_dict,data[t+1],t+1,N_ITER+1)
        if objective_value(case_max) < intradayMaxToken[t][:cost_to_go_copy] ||  N_ITER == 1
            om = intradayMax[t][:overestimator]
            v_upper = objective_value(case_max)
            @constraint(om,om[:μ] + sum(om[:τ_Pg][gen]*value(intraday[t][:Pg][gen]) for gen in keys(intradayMax[t][:Pg])) +
            sum(om[:τ_Ses][gen]*value(intraday[t][:Ses][gen]) for gen in keys(intradayMax[t][:Ses])) <=
            v_upper)#add cut to overestimator
        end
        intradayMax[t+1] = case_max_clean
        # **solve the updated lower problem**
        for (idx,gen) in enumerate(keys(case_dict[:windfarm]))
            fix(intraday[t+1][:Pw_err][gen],wst_vertex[idx]*α*data[t+1][:wind_power][gen])
        end
        optimize!(intraday[t+1])
        global total_solves
        total_solves += 1
        # **update the underestimator**
        lower_cut = AffExpr()
        for gen in keys(case_dict[:battery])
            #firstly the optimistic(base) case
            pi =  dual(FixRef(intraday[t+1][:Sesa][gen]))
            add_to_expression!(lower_cut,pi * (intraday[t][:Ses][gen] - value(intraday[t][:Ses][gen])))
            # pi_pgl = shadow_price()
            # add_to_expression!(dayahead_cut,pi_pgu*(dayahead[:Pgu][gen] - value(dayahead[:Pgu][gen])) + pi_pgl*(dayahead[:Pgu][gen] - value(dayahead[:Pgu][gen])))
        end
        for gen in intraday[t][:gens]
            pi =  dual(FixRef(intraday[t+1][:Pga][gen]))
            add_to_expression!(lower_cut,pi * (intraday[t][:Pg][gen] - value(intraday[t][:Pg][gen])))
        end
        # v_lower = sum(value(intraday[τ][:cost_now]) for τ in t+1:T)
        v_lower = objective_value(intraday[t+1])
        add_to_expression!(lower_cut,v_lower)
        @constraint(intraday[t],intraday[t][:cost_to_go] >= lower_cut)
    end
    return intraday,intradayMax
end

function stage_with_max_gap(intraday::Array{JuMP.Model},intradayMax::Array{JuMP.Model})
    # LsqFit.@. model(x, p) = p[1]*exp(-(x-1)*p[2])
    # xdata = 1:(length(intraday)-3)
    # ydata = [(objective_value(intradayMax[t]) - objective_value(intraday[t]))/objective_value(intradayMax[t]) for t in xdata]
    # p0 = [ydata[1],2]
    # fit = LsqFit.curve_fit(model, xdata, ydata, p0)
    # ydata_fitted = model(xdata,fit.param)
    # normalized = ydata./ydata_fitted
    # ~,start = findmin(normalized)
    # start = randperm(T)[1]
    return 1
end

function mature_stage(model::RealTimeDispatchModel,additional)
    if minimum(additional[:gaphourly][1:end-1]) < 0.01 && maximum(additional[:gaphourly][1:end-1]) < 0.2 
        stop = findfirst(x->x<=0.01,additional[:gaphourly][1:end-1])
    else 
        stop = length(model.intraday)
    end
    # @info("stop = $stop")
    return stop
end
function RDDP(model::RealTimeDispatchModel,config)
    # initialize remote solvers
    @info("$(length(workers())) workers really.")
    n_iter = 0
    global total_solves,sceanio_not_change_num
    total_solves = 0
    additional = Dict()
    solution_status = DataFrame(UpperBound=[],LowerBound=[],Gap=[],Time=[],TotalSolves=[]) 
    print_banner(stdout)
    print_iteration_header(stdout)
    start = 1
    stop = length(model.intraday)
    t1 = time()
    gap_temp = 9999
    no_improvement = 0
    sceanio_not_change_num = 0
    while true
        n_iter += 1
        additional = Problems.ForwardPassPrimal(model,n_iter,start,stop)
        if n_iter >= 5
            # stop = mature_stage(model,additional)
            # stop = max(stop,4)
        end
        # ______________________________________________________________________________________________________________
        
        additional[:Iteration] = n_iter
        additional[:TotalSolves] = total_solves
        t2 = time()
        additional[:Time] = t2 - t1
        print_iteration(stdout,additional)
        if additional[:Gap] >= gap_temp
            no_improvement += 1
        else
            no_improvement = 0
        end
        gap_temp = additional[:Gap]
        if n_iter >= 2
            push!(solution_status,additional)
        end
        if n_iter >= config["MaxIteration"] || additional[:Time] >= config["MaxTime"] || no_improvement >= config["no_improvement"] || additional[:Gap] < -0.01
            printstyled("Fail to converge. ";color=:red)
            print("$n_iter iterations in $(additional[:Time]) seconds. \n")
            break
        end
        if n_iter >=2 && (abs(additional[:Gap]) <= config["Gap"] || stop == 1)
            printstyled("Converged. ";color=:green)
            print("$n_iter iterations in $(round(additional[:Time];digits=3)) seconds. \n")
            break
        end
        @suppress_out Problems.BackwardPassPrimal(model,n_iter,start,stop)
    end
    return solution_status
end
function evalutaion(model::RealTimeDispatchModel,n=1000)
    # Random.seed!(1235)
    T = length(model.intraday)
    trajectory = []
    vts = model.vertice
    nk = length(vts[1])
    for i = 1:n
        trajectory_i = [vts[t][randperm(nk)[1]] for t in 1:T]
        push!(trajectory,trajectory_i)
    end
    objectives = []
    # p = Progress(n)
    for k in 1:n
        # total_cost[k] = value(dayahead[:cost_now])
        m = model.intraday
        total_cost = 0
        for t in 1:T
            fixall(model,t)
            for (idx,gen) in enumerate(keys(model.case_dict[:windfarm]))
                vertex = trajectory[k][t]
                windpower = model.data[t][:wind_power]
                @assert abs(vertex[idx])*m[t][:α] <= 0.5
                fix(m[t][:Pw_err][gen],vertex[idx]*m[t][:α]*windpower[gen])
            end
            @suppress_out optimize!(m[t])
            # println(value(m[t][:cost_now]))
            @assert termination_status(m[t]) == MOI.OPTIMAL #OPTIMAL
            total_cost += value(m[t][:cost_now])
        end
        push!(objectives,total_cost)
    end
    m = model.intraday
    worst = 0
    for t in 1:T
        fixall(model,t)
        for (idx,gen) in enumerate(keys(model.case_dict[:windfarm]))
            vertex = model.intradayMax[t][:wst_vertex]
            windpower = model.data[t][:wind_power]
            @assert abs(vertex[idx])*m[t][:α] <= 0.5
            fix(m[t][:Pw_err][gen],vertex[idx]*m[t][:α]*windpower[gen])
        end
        @suppress_out optimize!(m[t])
        # println(value(m[t][:cost_now]))
        @assert termination_status(m[t]) == MOI.OPTIMAL #OPTIMAL
        worst += value(m[t][:cost_now])
    end
    # while sum(objectives.>0) < n
    #     update!(p,sum(objectives.>0))
    # end
    return objectives,worst
end
function peaksOverThresholdEstimator(data::Array,quantile_α)
    u = quantile(data,quantile_α)
    extreme_values = filter(x->x>=u,data)
    nev = (extreme_values .- minimum(extreme_values))/std(extreme_values)
    N = length(nev)
    function f!(F,x)
        F[1] = 1/x[1] - (1/x[2] + 1) * 1/N * sum(nev[i]/(1 + x[1]*nev[i]) for i in 1:N)
        F[2] = 1/N * sum(log(1 + x[1]*nev[i]) for i in 1:N) - x[2]
    end
    x = NLsolve.nlsolve(f!,[0.1,0.1],autodiff=:forward)
    x = x.zero
    dist = Distributions.GeneralizedPareto(x[2]/x[1],x[2])
    return dist,minimum(extreme_values),std(extreme_values)
end
end

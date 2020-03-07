# using Revise
# using JuMP,CPLEX
# using JuMP,CDDLib,Polyhedra,Gurobi
# # includet("Problems.jl")

# # # optimize!(m,with_optimizer(Gurobi.Optimizer));
# # poly = polyhedron(m,CDDLib.Library());
# # # v = vrep(poly);
# # vertice = collect(points(poly))
# # # vertice = Problems.PolygonUncertaintySet([0,0,0],[1 0 0;0 1 0;0 0 1],1.5)
# # P = vrep([v for v in vertice])
# # Q = polyhedron(P)
# # removevredundancy!(Q)
# # simplex = hrep([HalfSpace([0,-1],-999),HalfSpace([-1,0],0),HalfSpace([1,0],100)])
# # poly = polyhedron(simplex,CDDLib.Library())
# # p = vrep([[50,50]])
# # simplex_updated = hrep(convexhull(poly,p))
# # poly = polyhedron(simplex_updated,CDDLib.Library())
# # p = vrep([[40,50]])
# # simplex_updated = hrep(convexhull(poly,p))
using Distributed,SharedArrays
addprocs(5)
@everywhere using JuMP, Gurobi
function make_model()
    m = Model(with_optimizer(Gurobi.Optimizer));
    @variable(m,x[1:1000]>=0)
    @constraint(m,sum(x)>=100 )
    @objective(m,Min,sum(x[i]^2*i for i in 1:1000))
    return m
end
# function solve_multicore(model)
#     Base.GC.enable(false)
#     Base.GC.enable(false)
#     envs = [Gurobi.Env() for i = 1:25]
#     envs_pointer = SharedArray{UInt64}(100)
#     for i = 1:25
#         envs_pointer[i] = UInt64(envs[i].ptr_env)
#     end
#     a = SharedArray{Float64}(100)
#     @sync @distributed for i = 1:100
#         p = Ptr{Nothing}(envs_pointer[myid()])
#         println(p)
#         env = envs[myid()]
#         env.ptr_env = p
#         set_optimizer(model[i],with_optimizer(Gurobi.Optimizer,env))
#         optimize!(model[i])
#         # # a[i] = objective_value(model[i])
#         # println(myid())
#     end
#     Base.GC.enable(true)
#     Base.GC.enable(true)
# end
m = make_model()
model = [m for i = 1:100]
# container = []
# # solve_multicore(model)
# for w in workers()
#     @spawnat w container
#     @spawnat w push!(container,Gurobi.Env())
# end
function f()
    @sync @distributed for i = 1:100
        println(Main.env)
        set_optimizer(model[i],with_optimizer(Gurobi.Optimizer,env))
        # optimize!(model[i])
        # # a[i] = objective_value(model[i])
        # println(myid())
    end
end    
@everywhere const env = Gurobi.Env()
a = SharedArray{Float64}(100)
least = minimum(workers())
f()
Base.GC.enable(true)
Base.GC.enable(true)

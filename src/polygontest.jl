# using Revise
# using JuMP,CPLEX
# using JuMP,CDDLib,Polyhedra,Gurobi
# # includet("Problems.jl")
# # m = Model();
# # @variable(m,0<=x[1:4]<=1);
# # @constraint(m,x[3]+x[4]>=+x[1]+x[2]);
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
using JuMP, Gurobi
model = Model(with_optimizer(Gurobi.Optimizer))
@variable(model, x)
@variable(model, y)
@variable(model, z)
@constraint(model, soccon, [x; y; z] in SecondOrderCone())
@constraint(model, eqcon, x == 1)
@objective(model, Min, y + z)
optimize!(model)
model2 = deepcopy(model)
@objective(model, Min, 2*y + z)
optimize!(model)

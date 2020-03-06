# syntax: proto3
using ProtoBuf
import ProtoBuf.meta

mutable struct OBJECTIVE <: ProtoType
    minFuelCost::Float64
    minCO2::Float64
    minERR::Float64
    minWindCurtailment::Float64
    OBJECTIVE(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct OBJECTIVE

mutable struct WIND <: ProtoType
    bus::Int32
    upper_bound::Base.Vector{Float64}
    lower_bound::Base.Vector{Float64}
    WIND(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct WIND
const __pack_WIND = Symbol[:upper_bound,:lower_bound]
meta(t::Type{WIND}) = meta(t, ProtoBuf.DEF_REQ, ProtoBuf.DEF_FNUM, ProtoBuf.DEF_VAL, true, __pack_WIND, ProtoBuf.DEF_WTYPES, ProtoBuf.DEF_ONEOFS, ProtoBuf.DEF_ONEOF_NAMES, ProtoBuf.DEF_FIELD_TYPES)

mutable struct DRFunction <: ProtoType
    c0::Float64
    c1::Float64
    c2::Float64
    DRFunction(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct DRFunction

mutable struct LOAD <: ProtoType
    bus::Int32
    _type::AbstractString
    load_of_a_day::Base.Vector{Float64}
    drfunction::DRFunction
    LOAD(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct LOAD
const __pack_LOAD = Symbol[:load_of_a_day]
meta(t::Type{LOAD}) = meta(t, ProtoBuf.DEF_REQ, ProtoBuf.DEF_FNUM, ProtoBuf.DEF_VAL, true, __pack_LOAD, ProtoBuf.DEF_WTYPES, ProtoBuf.DEF_ONEOFS, ProtoBuf.DEF_ONEOF_NAMES, ProtoBuf.DEF_FIELD_TYPES)

mutable struct HYDRO <: ProtoType
    bus::Int32
    initial_height::Base.Vector{Float64}
    HYDRO(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct HYDRO
const __pack_HYDRO = Symbol[:initial_height]
meta(t::Type{HYDRO}) = meta(t, ProtoBuf.DEF_REQ, ProtoBuf.DEF_FNUM, ProtoBuf.DEF_VAL, true, __pack_HYDRO, ProtoBuf.DEF_WTYPES, ProtoBuf.DEF_ONEOFS, ProtoBuf.DEF_ONEOF_NAMES, ProtoBuf.DEF_FIELD_TYPES)

mutable struct PSREQ <: ProtoType
    casefile::Array{UInt8,1}
    load_by_buses::Base.Vector{LOAD}
    wind_by_buses::Base.Vector{WIND}
    hydro_by_buses::Base.Vector{HYDRO}
    price::Base.Vector{Float64}
    objective::OBJECTIVE
    system_load::Base.Vector{Float64}
    PSREQ(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct PSREQ
const __fnum_PSREQ = Int[1,3,4,5,6,7,8]
const __pack_PSREQ = Symbol[:price,:system_load]
meta(t::Type{PSREQ}) = meta(t, ProtoBuf.DEF_REQ, __fnum_PSREQ, ProtoBuf.DEF_VAL, true, __pack_PSREQ, ProtoBuf.DEF_WTYPES, ProtoBuf.DEF_ONEOFS, ProtoBuf.DEF_ONEOF_NAMES, ProtoBuf.DEF_FIELD_TYPES)

mutable struct STATUS <: ProtoType
    gen::Int32
    status_of_a_day::Base.Vector{Int32}
    STATUS(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct STATUS
const __pack_STATUS = Symbol[:status_of_a_day]
meta(t::Type{STATUS}) = meta(t, ProtoBuf.DEF_REQ, ProtoBuf.DEF_FNUM, ProtoBuf.DEF_VAL, true, __pack_STATUS, ProtoBuf.DEF_WTYPES, ProtoBuf.DEF_ONEOFS, ProtoBuf.DEF_ONEOF_NAMES, ProtoBuf.DEF_FIELD_TYPES)

mutable struct POWER <: ProtoType
    gen::Int32
    power_of_a_day::Base.Vector{Float64}
    power_of_a_day_upper_bound::Base.Vector{Float64}
    power_of_a_day_lower_bound::Base.Vector{Float64}
    POWER(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct POWER
const __pack_POWER = Symbol[:power_of_a_day,:power_of_a_day_upper_bound,:power_of_a_day_lower_bound]
meta(t::Type{POWER}) = meta(t, ProtoBuf.DEF_REQ, ProtoBuf.DEF_FNUM, ProtoBuf.DEF_VAL, true, __pack_POWER, ProtoBuf.DEF_WTYPES, ProtoBuf.DEF_ONEOFS, ProtoBuf.DEF_ONEOF_NAMES, ProtoBuf.DEF_FIELD_TYPES)

mutable struct POWER_FLOW <: ProtoType
    f_bus::Int32
    t_bus::Int32
    power_flow_of_a_day::Base.Vector{Float64}
    POWER_FLOW(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct POWER_FLOW
const __pack_POWER_FLOW = Symbol[:power_flow_of_a_day]
meta(t::Type{POWER_FLOW}) = meta(t, ProtoBuf.DEF_REQ, ProtoBuf.DEF_FNUM, ProtoBuf.DEF_VAL, true, __pack_POWER_FLOW, ProtoBuf.DEF_WTYPES, ProtoBuf.DEF_ONEOFS, ProtoBuf.DEF_ONEOF_NAMES, ProtoBuf.DEF_FIELD_TYPES)

mutable struct PSRESP <: ProtoType
    termination_condition::AbstractString
    optimal_value::Float64
    status::Base.Vector{STATUS}
    power::Base.Vector{POWER}
    power_flow::Base.Vector{POWER_FLOW}
    error::AbstractString
    PSRESP(; kwargs...) = (o=new(); fillunset(o); isempty(kwargs) || ProtoBuf._protobuild(o, kwargs); o)
end #mutable struct PSRESP

# service methods for EconomicalRegion
const _EconomicalRegion_methods = MethodDescriptor[
        MethodDescriptor("GetEconomicalRegion", 1, PSREQ, PSRESP)
    ] # const _EconomicalRegion_methods
const _EconomicalRegion_desc = ServiceDescriptor("EconomicalRegion", 1, _EconomicalRegion_methods)

EconomicalRegion(impl::Module) = ProtoService(_EconomicalRegion_desc, impl)

mutable struct EconomicalRegionStub <: AbstractProtoServiceStub{false}
    impl::ProtoServiceStub
    EconomicalRegionStub(channel::ProtoRpcChannel) = new(ProtoServiceStub(_EconomicalRegion_desc, channel))
end # mutable struct EconomicalRegionStub

mutable struct EconomicalRegionBlockingStub <: AbstractProtoServiceStub{true}
    impl::ProtoServiceBlockingStub
    EconomicalRegionBlockingStub(channel::ProtoRpcChannel) = new(ProtoServiceBlockingStub(_EconomicalRegion_desc, channel))
end # mutable struct EconomicalRegionBlockingStub

GetEconomicalRegion(stub::EconomicalRegionStub, controller::ProtoRpcController, inp::PSREQ, done::Function) = call_method(stub.impl, _EconomicalRegion_methods[1], controller, inp, done)
GetEconomicalRegion(stub::EconomicalRegionBlockingStub, controller::ProtoRpcController, inp::PSREQ) = call_method(stub.impl, _EconomicalRegion_methods[1], controller, inp)

export PSREQ, OBJECTIVE, LOAD, WIND, DRFunction, HYDRO, PSRESP, STATUS, POWER, POWER_FLOW, EconomicalRegion, EconomicalRegionStub, EconomicalRegionBlockingStub, GetEconomicalRegion

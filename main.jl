using PowerModels
using DataStructures: SortedDict
using JuMP
using Gurobi
using DataFrames
using LinearAlgebra
using CSV
using Distributions

# load scripts
include("scripts/data_manager.jl")
include("scripts/fun_centralized_OPF.jl")
include("scripts/fun_compute_sensitivity.jl")
include("scripts/fun_consensus_update.jl")
include("scripts/fun_dual_update.jl")
include("scripts/fun_residual_update.jl")
include("scripts/fun_reveal_load.jl")
include("scripts/fun_voltage_update.jl")

# load data
caseID="testbeds/pglib_opf_case14_ieee.m"
(gen,bus,line,B,refbus)=load_data(caseID)
# initialize Gurobi environment
gurobi_env = Gurobi.Env()
# solve the centralized OPF problem
(cost_c,dispatch_c,power_flow_c)=OPF_centralized(gen,bus,line,B,refbus)
# create and specify ADMM and model parameters
ν̅ = 15000
μ = zeros(length(bus),length(bus),ν̅)
μ̃ = zeros(length(bus),length(bus),ν̅)
θ = zeros(length(bus),length(bus),ν̅)
θ̃ = zeros(length(bus),length(bus),ν̅)
θ̅ = zeros(length(bus),ν̅)
d = zeros(length(bus),length(bus),ν̅)
p = zeros(length(gen),ν̅)
l = zeros(length(bus),ν̅)
cost = zeros(1)
ρ = 1e3
γ = 1e-2
ν̃ = zeros(1)
# create and specify differential privacy parameters
ϵ = 1
α = 0.1
method = "PVP"
Δ_op=sensitivities(gen,bus,line,B,refbus,ρ,method,α)
ξ = zeros(length(bus),length(bus))
for i in 1:length(bus)
    for j in bus[i].N
        ξ[i,j,:] = rand(Laplace(0,Δ_op[i,j]/ϵ),1)
    end
end
# solve OPF using ADMM
for ν in 2:ν̅
    if method == "DVP"
        μ̃[:,:,ν] = μ[:,:,ν-1] .+ ξ[:,:]
        (θ[:,:,ν],cost[1],p[:,ν],l[:,ν]) = update_θ(gen,bus,line,B,refbus,μ̃[:,:,ν],θ̅[:,ν-1],ρ)
        d[:,:,ν] = reveal_load(bus,gen,B,ρ,μ[:,:,ν-1],θ̅[:,ν-1],θ[:,:,ν])
        θ̅[:,ν]=update_θ̅(bus,θ[:,:,ν])
        μ[:,:,ν]=update_μ(bus,ρ,θ[:,:,ν],θ̅[:,ν],μ[:,:,ν-1])
        Γ = residual(bus,θ[:,:,ν],θ̅[:,ν])
    end
    if method == "PVP"
        (θ[:,:,ν],cost[1],p[:,ν],l[:,ν]) = update_θ(gen,bus,line,B,refbus,μ[:,:,ν-1],θ̅[:,ν-1],ρ)
        θ̃[:,:,ν] = θ[:,:,ν] .+ ξ[:,:]
        d[:,:,ν] = reveal_load(bus,gen,B,ρ,μ[:,:,ν-1],θ̅[:,ν-1],θ̃[:,:,ν])
        θ̅[:,ν]=update_θ̅(bus,θ̃[:,:,ν])
        μ[:,:,ν]=update_μ(bus,ρ,θ̃[:,:,ν],θ̅[:,ν],μ[:,:,ν-1])
        Γ = residual(bus,θ̃[:,:,ν],θ̅[:,ν])
    end
    ν % 100 == 0 ? println("ν --- $(ν) ... res --- $(round(Γ,digits=5))") : NaN
    if Γ <= γ || ν == ν̅
        ν̃[1] = ν
        @info("ADMM terminates at iteration $(ν)")
        break
    end
end
# prepare results
load_inference = DataFrame(node=Any[],actual=Any[],observed=Any[])
for i in 1:length(bus)
    push!(load_inference,[i,bus[i].d,d[i,i,Int(ν̃[1])]])
end
node_dispatch = DataFrame(node=Any[],non_private=Any[],private=Any[])
for i in 1:length(bus)
    bus[i].type == 1 ? push!(node_dispatch,[i,dispatch_c[bus[i].G[1],3],p[bus[i].G[1],Int(ν̃[1])]]) : NaN
    bus[i].type == 2 ? push!(node_dispatch,[i,dispatch_c[i,4],l[i,Int(ν̃[1])]]) : NaN
end
flow_dispatch = DataFrame(line=Any[], b_f = Any[], b_t = Any[], flow_non_private = Any[], flow_private = Any[])
for l in 1:length(line)
    push!(flow_dispatch,[l,line[l].b_f,line[l].b_t,round(line[l].β * (dispatch_c[line[l].b_f,6]-dispatch_c[line[l].b_t,6]),digits=3),round(line[l].β * (θ̅[line[l].b_f,Int(ν̃[1])]-θ̅[line[l].b_t,Int(ν̃[1])]),digits=3)])
end
# print results
println("Optimality loss ---> $(abs(cost_c-cost[1])/cost_c*100)%")
println("Comparison of the non-private and differentially private load inference")
println(load_inference)
println("Comparison of the optimal and differentially private node dispatch")
println(node_dispatch)
println("Comparison of the optimal and differentially private flow dispatch")
println(flow_dispatch)

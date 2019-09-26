function sensitivities(gen,bus,line,B,refbus,ρ,method,α)
    """
    The function takes network data, method, and adjacency coefficient to compute
        the sensitivity of the agent voltage update to the value of loads. The
        method can be either "PVP" (primal variable perturbation) that uses
        formulation (7) in the paper or "DVP" (dual variable perturbation) that
        uses formulation (9) in the paper.
    """
    Ng = length(gen)
    Nb = length(bus)
    PV_buses=Set{Int64}()
    PQ_buses=Set{Int64}()
    for i in 1:Nb
        bus[i].type == 1 ? push!(PV_buses,i) : NaN
        bus[i].type == 2 ? push!(PQ_buses,i) : NaN
    end
    m = Model(with_optimizer(Gurobi.Optimizer,gurobi_env,Method=2,Presolve=1,OutputFlag=0))
    @variable(m, p[1:Ng]>=0)
    @variable(m, l[1:Nb])
    @variable(m, θ[1:Nb,1:Nb])
    @objective(m, Min, sum(gen[g].c2*p[g]^2 + gen[g].c1*p[g] + gen[g].c0 for g in 1:Ng) + sum(bus[i].c*l[i]^2 for i in 1:Nb) + ρ/2 * sum((-θ[i,j])^2 for i in 1:Nb for j in bus[i].N))
    @constraint(m, λ_PV[i=PV_buses], sum(B[i,j]*θ[i,j] for j in bus[i].N) == sum(p[g] for g in bus[i].G) - bus[i].d)
    @constraint(m, λ_PQ[i=PQ_buses], sum(B[i,j]*θ[i,j] for j in bus[i].N) == l[i] - bus[i].d)
    @constraint(m, κ[i=1:Nb], θ[i,refbus] == 0)
    optimize!(m)
    save_θ_α_0 = JuMP.value.(θ)
    for i=PV_buses
        delete(m, λ_PV[i])
    end
    for i=PQ_buses
        delete(m, λ_PQ[i])
    end
    @constraint(m, λ_PV_α[i=PV_buses], sum(B[i,j]*θ[i,j] for j in bus[i].N) == sum(p[g] for g in bus[i].G) - bus[i].d*(1-α))
    @constraint(m, λ_PQ_α[i=PQ_buses], sum(B[i,j]*θ[i,j] for j in bus[i].N) == l[i] - bus[i].d*(1-α))
    optimize!(m)
    save_θ_α = JuMP.value.(θ)
    δθ = zeros(Nb,Nb)
    for i in 1:Nb, j in bus[i].N
        if method == "PVP"
            δθ[i,j] = norm(save_θ_α_0[i,j]-save_θ_α[i,j]) + 1e-10
        end
        if method == "DVP" && bus[i].type == 1
            δθ[i,j] = norm(save_θ_α_0[i,j]-save_θ_α[i,j])*(2*gen[bus[i].G[1]].c2*B[i,j]^2+ρ) + 1e-10
        end
        if method == "DVP" && bus[i].type == 2
            δθ[i,j] = norm(save_θ_α_0[i,j]-save_θ_α[i,j])*(2*bus[i].c*B[i,j]^2+ρ) + 1e-10
        end
    end
    return δθ
end

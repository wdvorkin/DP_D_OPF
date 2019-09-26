function update_θ(gen,bus,line,B,refbus,μ,θ̅,ρ)
    """
    The function takes network data, dual, and consensus variables and returns
         local optimal power flow solution of each node according to the
         formulation (4a) in the paper.
    """
    Ng = length(gen)
    Nb = length(bus)
    Nl = length(line)
    model = Model(with_optimizer(Gurobi.Optimizer,gurobi_env,Method=2,Presolve=1,OutputFlag=0))
    @variable(model, p[1:Ng])
    @variable(model, θ[1:Nb,1:Nb])
    @variable(model, l[1:Nb])
    @constraint(model, ϕ̲[i=1:Nb,g=bus[i].G], gen[g].p̲ <= p[g])
    @constraint(model, ϕ̅[i=1:Nb,g=bus[i].G], p[g] <= gen[g].p̅)
    @constraint(model, ψ̲[i=1:Nb], -bus[i].d <= l[i])
    @constraint(model, ψ̅[i=1:Nb], l[i] <= bus[i].d)
    @constraint(model, η̲[i=1:Nb,l=bus[i].Λ], -line[l].f̅ <= line[l].β * (θ[i,line[l].b_f]-θ[i,line[l].b_t]))
    @constraint(model, η̅[i=1:Nb,l=bus[i].Λ],  line[l].β * (θ[i,line[l].b_f]-θ[i,line[l].b_t]) <= line[l].f̅)
    @constraint(model, λ[i=1:Nb], sum(B[i,j]*θ[i,j] for j in bus[i].N) == sum(p[g] for g in bus[i].G) - bus[i].d + l[i])
    @constraint(model, κ[i=1:Nb], θ[i,refbus] == 0)
    for i in 1:Nb
        bus[i].type == 1 ? @constraint(model, l[i] == 0) : NaN
    end
    @objective(model, Min, sum(gen[g].c2*p[g]^2 + gen[g].c1*p[g] + gen[g].c0 for g in 1:Ng) + sum(bus[i].c*l[i]^2 for i in 1:Nb) - sum(μ[i,j]*θ[i,j] for i in 1:Nb for j in bus[i].N) + ρ/2 * sum((θ̅[j]-θ[i,j])^2 for i in 1:Nb for j in bus[i].N))
    optimize!(model)
    cost = sum(gen[g].c2*JuMP.value(p[g])^2 + gen[g].c1*JuMP.value(p[g]) + gen[g].c0 for g in 1:Ng) + sum(bus[i].c*JuMP.value(l[i])^2 for i in 1:Nb)
    return JuMP.value.(θ),cost,JuMP.value.(p),JuMP.value.(l)
end

function OPF_centralized(gen,bus,line,B,refbus)
    """
    The function takes network data and computes optimal power flow
        solution using a centralized optimization according to the
        formulation (1) in the paper.
    """
    Ng = length(gen)
    Nb = length(bus)
    Nl = length(line)
    m = Model(with_optimizer(Gurobi.Optimizer, gurobi_env, Presolve=0, OutputFlag = 0))
    @variable(m, p[1:Ng])
    @variable(m, θ[1:Nb])
    @variable(m, l[1:Nb])
    # minimize cost
    @objective(m, Min, sum(gen[g].c2*p[g]^2 + gen[g].c1*p[g] + gen[g].c0 for g in 1:Ng) + sum(bus[i].c*l[i]^2 for i in 1:Nb))
    # DC-OPF constraints
    @constraint(m, ϕ[g=1:Ng], gen[g].p̲ <= p[g] <= gen[g].p̅)
    @constraint(m, ψ[i=1:Nb], 0 <= l[i] <= bus[i].d)
    @constraint(m, η[l=1:Nl], -line[l].f̅ <= line[l].β * (θ[line[l].b_f]-θ[line[l].b_t]) <= line[l].f̅)
    @constraint(m, λ[i=1:Nb], sum(B[i,j]*θ[j] for j in bus[i].N) == sum(p[g] for g in bus[i].G) - bus[i].d + l[i])
    @constraint(m, κ, θ[refbus] == 0)
    for i in 1:Nb
        bus[i].type == 1 ? @constraint(m, l[i] == 0) : NaN
    end
    optimize!(m)
    status = termination_status(m)
    @info("Centralized OPF terminates with status $(status)")
    cost=round(JuMP.objective_value(m),digits=5)
    dispatch = DataFrame(bus = Any[], d = Any[], gen = Any[], adj = Any[], adj_dual = Any[], ang = Any[])
    bus_gen = zeros(length(bus))
    for g in 1:Ng, b in 1:Nb
        gen[g].bus == b ? bus_gen[b] = JuMP.value(p[g]) : NaN
    end
    for i in 1:Nb
        push!(dispatch,[i,bus[i].d,round(bus_gen[i],digits=2),round(JuMP.value(l[i]),digits=2),round(JuMP.dual(ψ[i]),digits=3), round(JuMP.value(θ[i]),digits=3)])
    end
    power_flow = DataFrame(line=Any[], b_f = Any[], b_t = Any[], flow = Any[], cong = Any[])
    for l in 1:Nl
        abs(line[l].f̅-abs(line[l].β * (JuMP.value(θ[line[l].b_f])-JuMP.value(θ[line[l].b_t])))) <= 1e-2 ? status = 1 : status = 0
        push!(power_flow,[l,line[l].b_f,line[l].b_t,round(line[l].β * (JuMP.value(θ[line[l].b_f])-JuMP.value(θ[line[l].b_t])),digits=3),status])
    end
    return cost, dispatch, power_flow
end

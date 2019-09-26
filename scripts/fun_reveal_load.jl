function reveal_load(bus,gen,B,ρ,μ,θ̅,θ)
    """
    The function takes network data and agent communication signals and computes
        the vectors of estimated loads. Depending on the node type, it uses the
        formulation either (5a) or (5b) in the paper
    """
    Nb = length(bus)
    d = zeros(Nb,Nb)
    for i in 1:Nb
        for j in bus[i].N
            bus[i].type == 0 ? d[i,j] = 0 : NaN
            bus[i].type == 1 ? d[i,j] = (μ[i,j] + ρ*(θ̅[j]-θ[i,j]) - gen[bus[i].G[1]].c1*B[i,j])/(2*gen[bus[i].G[1]].c2*B[i,j]) - sum(B[i,k]*θ[i,k] for k in bus[i].N) : NaN
            bus[i].type == 2 ? d[i,j] = (μ[i,j] + ρ*(θ̅[j]-θ[i,j]))/(2*bus[i].c*B[i,j]) - sum(B[i,k]*θ[i,k] for k in bus[i].N) : NaN
        end
    end
    return d
end

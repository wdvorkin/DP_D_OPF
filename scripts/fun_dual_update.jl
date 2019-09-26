function update_μ(bus,ρ,θ,θ̅,μ)
    """
    The function takes all previous ADMM updates and computes the dual variables
        of the consensus constraint for each agent according to the formulation
        (4c) in the paper.
    """
    Nb = length(bus)
    mu = ones(Nb,Nb)
    for i in 1:Nb, j in bus[i].N
        mu[i,j] = μ[i,j] + ρ * (θ̅[j] - θ[i,j])
    end
    return mu
end

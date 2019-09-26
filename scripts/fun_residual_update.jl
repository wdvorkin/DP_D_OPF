function residual(bus,θ,θ̅)
    """
    The function updates the overall residual of the consensus constraint
        according to the formulation (4d) in the paper
    """
    Nb = length(bus)
    Γ = sum(norm(θ̅[j]-θ[i,j]) for i in 1:Nb for j in bus[i].N)
    return Γ
end

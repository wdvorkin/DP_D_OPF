function update_θ̅(bus,θ)
    """
    The function updates consensus variables based on the agent voltage angle
        estimates across the neighborhood of each node. The computation of
        consensus variables is reduced to simple averaging.
        For details, see Section 7.1 in
            Boyd, Stephen, et al. "Distributed optimization and statistical
            learning via the alternating direction method of multipliers."
            Foundations and Trends® in Machine learning 3.1 (2011): 1-122.
    """
    Nb = length(bus)
    θ̅ = zeros(Nb)
    for i in 1:Nb
        θ̅[i] = sum(θ[j,i] for j in bus[i].N)/length(bus[i].N)
    end
    return θ̅
end

"""
The script builds the function to return network data of chosen test case
    from PowerModels.jl
"""
mutable struct  Bus
   ind::Int
   d::Float64
   c::Float64
   type::Int
   N::Vector{Int}
   G::Vector{Int}
   Λ::Vector{Int}
   function Bus(ind,d,c,type,N,G,Λ)
      i = new()
      i.ind  = ind
      i.d = d
      i.c = c
      i.type = type
      i.N = N
      i.G = G
      i.Λ = Λ
      return i
   end
end
mutable struct  Generator
   ind::Int
   c2::Any
   c1::Any
   c0::Any
   p̅::Any
   p̲::Any
   bus::Int
   function Generator(ind,c2,c1,c0,p̅,p̲,bus)
      i = new()
      i.ind  = ind
      i.c2 = c2
      i.c1 = c1
      i.c0 = c0
      i.p̅ = p̅
      i.p̲ = p̲
      i.bus = bus
      return i
   end
end
mutable struct  Line
   ind::Int
   b_f::Int
   b_t::Int
   β::Float64
   f̅::Float64
   θ̅::Float64
   θ̲::Float64
   function Line(ind,b_f,b_t,β,f̅,θ̅,θ̲)
      i = new()
      i.ind  = ind
      i.b_f = b_f
      i.b_t = b_t
      i.β = β
      i.f̅ = f̅
      i.θ̅ = θ̅
      i.θ̲ = θ̲
      return i
   end
end

function load_data(caseID)
    PowerModels.silence()
    data = PowerModels.parse_file(caseID)

    gen = Dict()
    for i in 1:length(data["gen"]) # remove generators with zero active power capacity
        data["gen"][string(i)]["pmax"]*data["baseMVA"] == 0 ? delete!(data["gen"],string(i)) : NaN
    end
    gen_key=collect(keys(data["gen"]))
    savec2 = zeros(length(gen_key))
    for i in 1:length(gen_key)
        ind = i
        if sum(data["gen"][gen_key[i]]["cost"]) > 0
            data["gen"][gen_key[i]]["ncost"] == 3 ? c2 = data["gen"][gen_key[i]]["cost"][1]/data["baseMVA"]^2 : c2 = 0
            data["gen"][gen_key[i]]["ncost"] == 3 ? c1 = data["gen"][gen_key[i]]["cost"][2]/data["baseMVA"] : c1 = data["gen"][gen_key[i]]["cost"][1]/data["baseMVA"]
            data["gen"][gen_key[i]]["ncost"] == 3 ? c0 = data["gen"][gen_key[i]]["cost"][3]  : c0 = 0
        else
            c2 = 0
            c1 = 0
            c0 = 0
        end
        c2 == 0 ? c2 = 0.01 * c1 : NaN
        p̅ = data["gen"][gen_key[i]]["pmax"]*data["baseMVA"]
        p̲ = data["gen"][gen_key[i]]["pmin"]*data["baseMVA"]
        bus = data["gen"][gen_key[i]]["gen_bus"]
        add_gen = Generator(ind,c2,c1,c0,p̅,p̲,bus)
        gen[add_gen.ind] = add_gen
        savec2[i] = c2
    end
    c2max = maximum(savec2)
    gen=SortedDict(gen)

    line = Dict()
    for l in 1:length(data["branch"])
        ind = l
        b_f = data["branch"][string(l)]["f_bus"]
        b_t = data["branch"][string(l)]["t_bus"]
        β = -imag(1/(data["branch"][string(l)]["br_r"] + data["branch"][string(l)]["br_x"]im))
        f̅ = data["branch"][string(l)]["rate_a"]*data["baseMVA"]
        θ̅ = data["branch"][string(l)]["angmax"]
        θ̲ = data["branch"][string(l)]["angmin"]
        add_line = Line(ind,b_f,b_t,β,f̅,θ̅,θ̲)
        line[add_line.ind] = add_line
    end
    line=SortedDict(line)

    bus = Dict()
    for b in 1:length(data["bus"])
        ind = b
        d = 0
        c = 10
        # *c2max
        type = 0
        for l in 1:length(data["load"])
            data["load"][string(l)]["load_bus"] == b ? d += data["load"][string(l)]["pd"]*data["baseMVA"] : NaN
        end
        N = Int[]
        G = Int[]
        Λ = Int[]
        add_bus = Bus(ind,d,c,type,N,G,Λ)
        bus[add_bus.ind] = add_bus
    end
    for g in 1:length(data["gen"]), b in 1:length(data["bus"])
        data["gen"][gen_key[g]]["gen_bus"] == b ? push!(bus[b].G,g) : NaN
        data["gen"][gen_key[g]]["gen_bus"] == b && data["gen"][gen_key[g]]["pmax"] > 0 ? bus[b].type = 1 : NaN
    end
    for b in 1:length(data["bus"])
        if bus[b].type != 1
            bus[b].d != 0 ? bus[b].type = 2 : NaN
        end
    end
    for l in 1:length(data["branch"]), b in 1:length(data["bus"])
        data["branch"][string(l)]["f_bus"] == b ? push!(bus[b].Λ,l) : NaN
        data["branch"][string(l)]["t_bus"] == b ? push!(bus[b].Λ,l) : NaN
    end
    for b in 1:length(data["bus"]), l in 1:length(data["branch"])
        data["branch"][string(l)]["f_bus"] == b ? push!(bus[b].N,data["branch"][string(l)]["t_bus"]) : NaN
        data["branch"][string(l)]["t_bus"] == b ? push!(bus[b].N,data["branch"][string(l)]["f_bus"]) : NaN
        push!(bus[b].N,b)
        bus[b].N = unique(bus[b].N)
    end
    bus=SortedDict(bus)

    Z_eq = zeros(Complex{Float64},length(bus),length(bus))
    B = zeros(length(bus),length(bus))
    for l in 1:length(line)
        Z_eq[line[l].b_f,line[l].b_t] == 0 ? Z_eq[line[l].b_f,line[l].b_t] = data["branch"][string(l)]["br_r"] + data["branch"][string(l)]["br_x"]im : Z_eq[line[l].b_f,line[l].b_t] = (Z_eq[line[l].b_f,line[l].b_t]*(data["branch"][string(l)]["br_r"] + data["branch"][string(l)]["br_x"]im))/(Z_eq[line[l].b_f,line[l].b_t]+(data["branch"][string(l)]["br_r"] + data["branch"][string(l)]["br_x"]im))
        Z_eq[line[l].b_t,line[l].b_f] == 0 ? Z_eq[line[l].b_t,line[l].b_f] = data["branch"][string(l)]["br_r"] + data["branch"][string(l)]["br_x"]im : Z_eq[line[l].b_t,line[l].b_f] = (Z_eq[line[l].b_t,line[l].b_f]*(data["branch"][string(l)]["br_r"] + data["branch"][string(l)]["br_x"]im))/(Z_eq[line[l].b_t,line[l].b_f]+(data["branch"][string(l)]["br_r"] + data["branch"][string(l)]["br_x"]im))
    end
    for i in 1:length(bus), j in 1:length(bus)
        Z_eq[i,j] != 0 ? B[i,j] = imag(1/Z_eq[i,j]) : 0
    end
    for i in 1:length(bus)
        a = -sum(B[i,:])
        B[i,i] = a
    end
    refbus = 1
    return gen,bus,line,B,refbus
end

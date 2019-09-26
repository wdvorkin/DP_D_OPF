# E-Companion for "Differentially Private Distributed Optimal Power Flow"
by  V. Dvorkin, P. Van Hentenryck, J. Kazempour and P. Pinson. 

The repository contains data and codes for the centralized OPF in (1) and two privacy-preserving ADMM Algorithms 1 and 2. The code returns the comparison of the centralized OPF solution with the solution of one of the privacy-preserving algorithms. The optimization models have been implemented by using [JuMP](https://github.com/JuliaOpt/JuMP.jl) in the [Julia](http://julialang.org/downloads/) programming language (v1.1.1.1). The following packages are required:
- [DataFrames.jl](https://github.com/DataFrames.jl/stable/)
- [DataStructures.jl](https://github.com/JuliaCollections/DataStructures.jl)
- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)
- [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)
- [JuMP](https://github.com/JuliaOpt/JuMP.jl)
- [LinearAlgebra.jl](https://github.com/JuliaStdlibs/LinearAlgebra.jl)
- [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl)

To run the program, make use of particular versions of these Julia packages with
```
julia> Pkg.pin("DataFrames", v"0.19.2")
julia> Pkg.pin("DataStructures", v"0.17.0")
julia> Pkg.pin("Distributions", v"0.21.1")
julia> Pkg.pin("Gurobi", v"0.6.0")
julia> Pkg.pin("JuMP", v"0.19.2")
julia> Pkg.pin("PowerModels", v"0.12.2")
```

mutable struct CriticalRegion 
    AS::Vector{Int16}
    Ath::Matrix{Float64}
    bth::Vector{Float64}
    x::Matrix{Float64}
    lam::Matrix{Float64}
    th::Vector{Float64}
end

struct MPLDP 
    MM::Matrix{Float64}
    M::Matrix{Float64}
    MRt::Matrix{Float64}
    RinvV::Matrix{Float64}
    d::Matrix{Float64}
    n_theta::Int64
    n::Int64
    bounds_table::Vector{Int64}
    norm_factors::Vector{Float64}
end

Base.@kwdef mutable struct Settings 
    eps_zero::Float64 = 1e-12
    verbose::Int64 = 1 
    store_AS::Bool = true
    store_points::Bool = false
    store_regions::Bool = true
    store_dual::Bool = false
    remove_redundant::Bool = true
    time_limit::Int64 = 1e5
    chunk_size::Int64 = 1e3 
    factorization::Symbol = :chol
    postcheck_rank::Bool = true
    lowdim_tol::Float64 = 0 
end

Settings(opts::Nothing) = Settings()
Settings(opts::Settings) = opts
function Settings(opts::Dict)
    out = Settings()
    for (key,value) in opts
        if hasproperty(out, Symbol(key))
            setproperty!(out, Symbol(key), value)
        else
            @warn "$key is not a valid setting"
        end
    end
    return out
end

mutable struct Workspace{T<:Integer}
    Ath::Matrix{Float64}
    bth::Vector{Float64}
    bth_lower::Vector{Float64}
    sense::Vector{Cint}
    m::Int64
    m0::Int64 
    DAQP_workspace::Ptr{Cvoid}
    ASs::BitMatrix
    nLPs::Int64
    S::Vector{T}
    Sdown::Vector{T}
    Sup::Vector{T}
    F::Vector{CriticalRegion}
    explored::Set{T}
    opts::Settings
    IS::BitVector
    AS::BitVector
    nAS::Int64
    norm_factors:: Vector{Float64}
end

mutable struct CriticalRegion 
    AS::Vector{Int16}
    Ath::Matrix{Float64}
    bth::Vector{Float64}
    xTH::Matrix{Float64}
    xC::Vector{Float64}
    th::Vector{Float64}
end

struct MPLDP 
    MM::Matrix{Float64}
    M::Matrix{Float64}
    MRt::Matrix{Float64}
    vTH::Matrix{Float64}
    vC::Vector{Float64}
    d::Matrix{Float64}
    n_theta::Int64
    n::Int64
    bounds_table::Vector{Int64}
end

Base.@kwdef mutable struct EMPCSettings 
    eps_zero::Float64 = 1e-12
    eps_gap::Float64  = 1e-6-1e-12
    verbose::Int64 = 1 
    store_AS::Bool = true
    store_points::Bool = false
    store_regions::Bool = true
    remove_redundant::Bool = true
    time_limit::Int64 = 1e5
    chunk_size::Int64 = 1e3 
    factorization::Symbol = :chol
    postcheck_rank::Bool = true
    lowdim_tol::Float64 = 0 
end

mutable struct EMPCWorkspace{T<:Integer}
    Ath::Matrix{Float64}
    bth::Vector{Float64}
    bth_lower::Vector{Float64}
    sense::Vector{Cint}
    m::Int64
    m0::Int64 
    DAQP_workspace::Ptr{Cvoid}
    ASs::BitMatrix
    nLPs::Int64
    S::Vector{Tuple{T,Int32}}
    Sdown::Vector{Tuple{T,Int32}}
    Sup::Vector{Tuple{T,Int32}}
    F::Vector{CriticalRegion}
    explored::Set{T}
    opts::EMPCSettings
    IS::BitVector
    AS::BitVector
    nAS::Int64
    norm_factors:: Vector{Float64}
    parents::Vector{Int32}
end

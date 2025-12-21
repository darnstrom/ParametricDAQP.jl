mutable struct CriticalRegion 
    AS::Vector{Int16}
    Ath::Matrix{Float64}
    bth::Vector{Float64}
    z::Matrix{Float64}
    lam::Matrix{Float64}
    th::Vector{Float64}
end

struct MPLDP 
    AHinvA::Matrix{Float64}
    M::Matrix{Float64}
    HinvAt::Matrix{Float64}
    HinvF::Matrix{Float64}
    d::Matrix{Float64}
    n_theta::Int64
    n::Int64
    bounds_table::Vector{Int64}
    norm_factors::Vector{Float64}
    eq_ids::Vector{Int64}
    out_lims::Matrix{Float64}
end
struct MPQP
    H::AbstractMatrix
    F::Matrix{Float64}
    A::Matrix{Float64}
    B::Matrix{Float64}
    n_theta::Int64
    n::Int64
    bounds_table::Vector{Int64}
    norm_factors::Vector{Float64}
    eq_ids::Vector{Int64}
    rank_defficiency::Int64
    out_inds::Vector{Int64}
    out_lims::Matrix{Float64}
end
struct MPVI
    # VI(Hx + Fθ + f, Ax <= Bθ + b)
    # where f, b are the last rows of F,B, respect.
    H::AbstractMatrix # size = n_x * n_x
    F::Matrix{Float64} # transpose: size = n_θ + 1 * n_x 
    A::Matrix{Float64} # size = n_constr * n_x
    B::Matrix{Float64} # transpose: size = n_θ + 1 * n_constr 
    AHinv::Matrix{Float64} # size = n_constr * n_x
    AHinvA::Matrix{Float64} # size = n_constr * n_constr
    n_theta::Int64
    n::Int64
    bounds_table::Vector{Int64}
    norm_factors::Vector{Float64}
    eq_ids::Vector{Int64} # introduced for compatibility reasons - only inequalities are considered in the VI case
    out_inds::Vector{Int64}
    out_lims::Matrix{Float64}
end

"""
    MPVI(
        H::AbstractMatrix{Float64},
        G::AbstractMatrix{Float64},
        f::AbstractVector{Float64},
        A::AbstractMatrix{Float64},
        E::AbstractMatrix{Float64},
        b::AbstractVector{Float64}
    )

Constructs an `MPVI` (Multi-Parametric Variational Inequality) problem instance.

# Arguments
- `H::AbstractMatrix{Float64}`: Matrix representing the linear part of the affine mapping.
- `G::AbstractMatrix{Float64}`: Matrix for the parametric part of the affine mapping.
- `f::AbstractVector{Float64}`: Vector for the constant part of the affine mapping.
- `A::AbstractMatrix{Float64}`: Constraint matrix.
- `E::AbstractMatrix{Float64}`: Matrix for the parametric part of the constraint right-hand side.
- `b::AbstractVector{Float64}`: Vector for the constant part of the constraint right-hand side.

# Returns
- An `MPVI`  if H is nonsymmetric
- An 'MPQP'  if H is symmetric and singular
- An 'MPLDP' if H is symmetric and nonsingular 

"""
function MPVI(
        H::AbstractMatrix{Float64},
        G::AbstractMatrix{Float64},
        f::AbstractVector{Float64},
        A::AbstractMatrix{Float64},
        E::AbstractMatrix{Float64},
        b::AbstractVector{Float64}
    )
    return setup_mpp((H=H,F=G,f=f,A=A,B=E,b=b))
end

Base.@kwdef mutable struct Settings 
    eps_zero::Float64 = 1e-12
    verbose::Int64 = 1 
    store_AS::Bool = true
    store_points::Bool = true
    store_regions::Bool = true
    store_dual::Bool = false
    remove_redundant::Bool = true
    time_limit::Int64 = 1e5
    region_limit::Int64 = 1e12
    chunk_size::Int64 = 1e3 
    factorization::Symbol = :chol
    postcheck_rank::Bool = true
    lowdim_tol::Float64 = 1e-12
    daqp_settings = Dict{Symbol,Any}()
end

Settings(opts::Nothing) = Settings()
Settings(opts::Settings) = opts
function Settings(opts::AbstractDict)
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
    DAQP_workspace::Ptr{DAQPBase.Workspace}
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
    x::Matrix{Float64}
end

struct Solution
    problem::Union{MPLDP,MPQP,MPVI}
    CRs::Vector{CriticalRegion}
    scaling::Vector{Float64} 
    translation::Vector{Float64}
    settings::Settings
    status::Symbol
end

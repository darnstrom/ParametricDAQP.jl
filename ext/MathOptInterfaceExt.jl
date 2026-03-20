module MathOptInterfaceExt

using ParametricDAQP
using MathOptInterface
using LinearAlgebra

include("./moi_bridge.jl")

function pdaqp_get_varids(idxmap::MathOptInterface.Utilities.IndexMap,vars)
    return [idxmap[var.index].value for var in vars]
end
function moi2mpqp(moi_model::MOI.ModelLike, vars; out_vars=nothing, eliminate_equalities=false)
    opt = MathOptInterface.instantiate(ParametricOptimizer;with_cache_type=Float64, with_bridge_type=Float64)
    MathOptInterface.copy_to(opt,moi_model)
    src, dest = opt.model, opt.model.optimizer
    idxmap = MathOptInterface.Utilities.IndexMap(dest, src) 
    check_attributes(dest,src)
    assign_constraint_rows!(dest, idxmap, src)

    A, bupper, blower, sense = process_constraints(dest, src, idxmap)
    H, f, c = process_objective(dest,src, idxmap)

    n,m = length(f),length(bupper)
    par_ids = pdaqp_get_varids(idxmap,vars)
    dec_ids = setdiff(1:n,par_ids)
    out_inds = isnothing(out_vars) ? dec_ids :  pdaqp_get_varids(idxmap,out_vars)
    out_inds = [findfirst(==(i),dec_ids) for i in out_inds]

    Htot,Ftot,ftot = H[dec_ids,dec_ids],H[dec_ids,par_ids],f[dec_ids]
    H,F,f = Htot,Ftot,ftot

    th_lb,th_ub = blower[par_ids],bupper[par_ids]
    th_lb[th_lb .< -1e20] .= -100
    th_ub[th_ub .> 1e20] .=  100

    mask = falses(m)
    mask[dec_ids] .= true
    mask[n+1:end] .= true
    mask .&=  (sense .!=DAQPBase.IMMUTABLE)

    B = Matrix([zeros(length(dec_ids),length(par_ids)); -A[par_ids,:]'])
    A = Matrix(A[dec_ids,:]')
    bl, bu = blower[mask],bupper[mask] 
    senses = sense[mask]

    mpp = (H=H,F=F,f=f,A=A,B=B,bl=bl,bu=bu,senses=senses)

    if eliminate_equalities
        mpp = ParametricDAQP.preprocess_eliminate_equalities(mpp)
    end
    TH = (lb=th_lb, ub = th_ub)
    return mpp,TH
end

function ParametricDAQP.mpsolve(model::MOI.ModelLike,vars;
        out_vars=nothing,opts=nothing, eliminate_equalities=true)
    mpqp,TH = moi2mpqp(model,vars;out_vars,eliminate_equalities)
    mpqp = ParametricDAQP.make_singlesided(mpqp)
    return ParametricDAQP.mpsolve(mpqp,TH;opts)
end

function ParametricDAQP.get_mpp(moi_model::MOI.ModelLike, vars; out_vars=nothing, eliminate_equalities=false)
    return moi2mpqp(moi_model,vars;out_vars,eliminate_equalities)
end

end

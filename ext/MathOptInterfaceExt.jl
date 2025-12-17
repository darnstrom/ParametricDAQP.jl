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
    #H+=1e-6I


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


    Aext = [I(n); A'][mask,:]
    bu = bupper[mask] 
    bl = blower[mask] 
    senses = sense[mask]

    mask_eq = senses .== DAQPBase.EQUALITY
    Aeq = Aext[mask_eq,:]
    beq = bu[mask_eq,:] 

    if(eliminate_equalities)
        # TODO: retrieve solution (Maybe this could be moved...)
        Weq = -Aeq[:,par_ids]
        Aeq = Aeq[:,dec_ids]
        z0,zTH,Z = Aeq\beq, Aeq\Weq, nullspace(Aeq)
        f,F = Z'*(f+H*z0), Z'*(F+H*zTH) 
        H = Z'H*Z
        Aineq = Aext[.~mask_eq,:]
        A,W = Aineq[:,dec_ids], -Aineq[:,par_ids]
        blineq,buineq =  bl[.~mask_eq] - A*z0, bu[.~mask_eq] - A*z0
        W, Aineq = W-A*zTH, A*Z

        Atot = [Aineq;-Aineq]
        btot = [buineq;-blineq]
        Wtot = [W;-W]
        mi = length(buineq)
        bounds_table = [mi+1:2*mi;1:mi]
        senses = repeat(senses[.~mask_eq],2)

    else
        Aineq = kron([1;-1],Aext[.~mask_eq,:])
        bineq = [bu[.~mask_eq];-bl[.~mask_eq]]

        mi,me = Int(length(bineq)/2),length(beq)
        bounds_table = [1:me; me+mi+1:me+2*mi;me+1:me+mi]
        senses = [senses[mask_eq];repeat(senses[.~mask_eq],2)]

        Atot = [Aeq;Aineq]
        btot = [beq;bineq]

        Wtot = -Atot[:,par_ids] 
        Atot = Atot[:,dec_ids] 
    end

    # Prune superluous constraints
    rm_ids = findall(btot[:] .>= 1e20)
    if(!isempty(rm_ids))
        bounds_table[bounds_table[rm_ids]] = bounds_table[rm_ids] # Make other bound point to itself
        #Correct bounds table 
        rm_offset, keep_ids = 1, Int[]
        for i in 1:length(btot)
            if(rm_offset <= length(rm_ids) && i==rm_ids[rm_offset])
                rm_offset+=1
            else
                bounds_table[i] -= (rm_offset-1)
                push!(keep_ids,i)
            end
        end
        Atot,btot,Wtot = Atot[keep_ids,:],btot[keep_ids],Wtot[keep_ids,:]
        senses,bounds_table = senses[keep_ids],bounds_table[keep_ids]
    end

    mpQP = (H=H,F=F,f=f,
            A = Atot, b=btot, W=Wtot, senses = senses, 
            bounds_table = bounds_table, out_inds = out_inds)
    TH = (lb=th_lb, ub = th_ub)
    return mpQP,TH
end


function ParametricDAQP.mpsolve(model::MOI.ModelLike,vars;out_vars=nothing,opts=nothing)
    mpqp,TH = moi2mpqp(model,vars;out_vars)
    return ParametricDAQP.mpsolve(mpqp,TH;opts)
end

end

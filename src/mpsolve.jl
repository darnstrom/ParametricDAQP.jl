## mpsolve 
# Compute the explicit solution to multi-parameteric QP
function mpsolve(mpQP,Θ;opts=nothing, AS0 = nothing) # bounds_table as option
    opts = Settings(opts)
    # Hande bounds 
    Δb = Θ.ub-Θ.lb 
    if any(Δb .< -opts.eps_zero) 
        opts.verbose > 0 && @error "Θ is empty"
        return nothing,nothing
    end


    # Check if some parameters are fixed
    fix_ids = findall(Δb .< opts.eps_zero) # ub[fix_ids] = lb[fix_ids]
    fix_vals = Θ.ub[fix_ids]
    if(!isempty(fix_ids))
        opts.verbose > 0 && @warn "θ$fix_ids fixed at $fix_vals"
        free_ids = setdiff(1:length(Θ.ub),fix_ids)
        if hasproperty(Θ,:A)
            Θ=(A=Θ.A[free_ids,:], b = Θ.b+Θ.A[fix_ids,:]'*fix_vals, 
               lb=Θ.lb[free_ids],ub=Θ.ub[free_ids])
        else
            Θ=(lb=Θ.lb[free_ids],ub=Θ.ub[free_ids])
        end
    end

    # Setup parametric problem and normalize
    prob = setup_mpp(mpQP;fix_ids,fix_vals)
    prob, Θ, tf = normalize_parameters(prob,Θ)
    if(isnothing(prob)) # Θ makes the problem trivially infeasible
        @warn "The parameter region of interest is empty"
        F,info = CriticalRegion[], (solve_time = 0, nCR = 0, nLPs = 0, nExplored = 0,status=:EmptyParameterRegion)
    end


    if(isnothing(AS0))
        AS0 = compute_AS0(prob,Θ)
    end
    if(isnothing(AS0))
        F,info = CriticalRegion[], (solve_time = 0, nCR = 0, nLPs = 0, nExplored = 0,status=:NoFeasibleParameter)
    else
        F,info = mpdaqp_explicit(prob,Θ,AS0;opts)
    end
    return Solution(F,tf.scaling,tf.center,opts,info.status), info
end
## Method based on combinatorial adjacency 
function mpdaqp_explicit(prob,Θ,AS0;opts = Settings())
    time_limit = opts.time_limit*1e9;
    status = :Solved
    t0  = time_ns()
    m = length(prob.norm_factors);

    # Handle zero rows
    id_cands = findall(prob.norm_factors .> opts.eps_zero)
    if(length(id_cands) < m)
        id_zeros = findall(prob.norm_factors .≤ opts.eps_zero)
        opts.verbose >  0 && @warn "Rows $id_zeros in A are zero → seen as parameter constraints" 
        if hasproperty(Θ,:A)
            A = [Θ.A -prob.d[1:end-1,id_zeros]]
            b = [Θ.b; prob.d[end,id_zeros]]
        else
            A = -prob.d[1:end-1,id_zeros]
            b = prob.d[end,id_zeros]
        end
        Θ = (A= A, b = b,lb = Θ.lb, ub = Θ.ub) 
    end

    setdiff!(id_cands,prob.eq_ids) # equality constraints cannot be modified

    # Initialize
    ws=setup_workspace(Θ,m;opts);
    as0 = as2uint(AS0,eltype(ws.S))
    push!(ws.S,as0)
    push!(ws.explored,as0)
    j = 0


    # Start exploration 
    while(!isempty(ws.S) || !isempty(ws.Sdown) || !isempty(ws.Sup))

        # Check time limit 
        if(time_ns()-t0 > time_limit)
            status = :TimeLimitReached
            break
        end

        while(length(ws.S) < opts.chunk_size && !isempty(ws.Sdown)) # First try to move down...
            explore_supersets(pop!(ws.Sdown),ws,id_cands,ws.S,prob.bounds_table)
        end
        while(length(ws.S) < opts.chunk_size && !isempty(ws.Sup)) # ... then try to move up
            explore_subsets(pop!(ws.Sup),ws,id_cands,ws.S)
        end

        opts.verbose>0 && print_ws(ws,(j+=1))

        # Process pending AS
        while(!isempty(ws.S))
            as = pop!(ws.S);
            region,up,down= isoptimal(as,ws,prob,opts)
            !isnothing(region) && push!(ws.F,region)
            up && push!(ws.Sup,as)
            down && push!(ws.Sdown,as)
        end

    end
    # Exploration completed, cleanup 
    DAQP.free_c_workspace(ws.DAQP_workspace)
    solve_time = (time_ns()-t0)/1e9
    opts.verbose>0 && print_final(ws)
    return ws.F, (solve_time = solve_time, nCR = length(ws.F), 
                  nLPs = ws.nLPs, nExplored = length(ws.explored),
                  status=status)
end
## Check if optimal
function isoptimal(as,ws,prob,opts)
    # Extract AS and IS for current candidate
    uint2as(as,ws,prob.bounds_table)
    ws.nAS > prob.n  && return nothing,true,false # LICQ trivially broken

    # Compute λ and μ 
    up, down = compute_λ_and_μ(ws,prob,opts)
    (up || down) && return nothing,up,down # Degeneracy

    # Reset feasibility workspace
    reset_workspace(ws) 
    # Solve primal+dual feasibility problem 
    normalize_model(ws;eps_zero=opts.eps_zero) || return nothing,false,false;

    ws.nLPs+=1
    if isfeasible(ws.DAQP_workspace, ws.m, 0)
        return extract_CR(ws,prob),true,ws.nAS ≤ prob.n
    else
        return nothing,false,false
    end
end
## Update optimization model
function normalize_model(ws;eps_zero=1e-12)
    norm_factor=0.0;
    # add λ ≥ 0 to feasibility model
    for i in (1+ws.m0: ws.m0+length(ws.AS))
        norm_factor = norm(view(ws.Ath,:,i),2);
        ws.norm_factors[i-ws.m0] = norm_factor
        if(norm_factor>eps_zero)
            rdiv!(view(ws.Ath,:,i),norm_factor)
            ws.bth[i]/=norm_factor
            ws.sense[i] = 0
        else
            (ws.bth[i]<-eps_zero) && return false # trivially infeasible 
            ws.sense[i] = 4 
        end
    end
    ws.m=length(ws.AS)+ws.m0;
    return true
end
## Compute slacks 
function compute_λ_and_μ(ws,prob::MPLDP,opts)
    if(opts.factorization == :qr)
        R = (ws.nAS > 0) ? UpperTriangular(qr(prob.M[ws.AS,:]').R) : UpperTriangular(zeros(0,0))
        if(any(abs(R[i,i]) <1e-12 for i in 1:ws.nAS))
            return true,false
        end
    else
        C = cholesky(prob.MM[ws.AS,ws.AS],check=false)
        if(!issuccess(C))
            return true, false
        end
        R = UpperTriangular(C.factors)
    end

    # Compute λ
    λTH = @view ws.Ath[:,ws.m0+1:ws.m0+ws.nAS]
    λC = @view ws.bth[ws.m0+1:ws.m0+ws.nAS]
    λTH .= @view prob.d[1:end-1,ws.AS]; rdiv!(rdiv!(λTH, R), adjoint(R))
    λC .= -@view prob.d[end,ws.AS]; ldiv!(R, ldiv!(adjoint(R), λC))
    # Compute μ
    μTH = @view ws.Ath[:,ws.m0+ws.nAS+1:end]
    μC = @view ws.bth[ws.m0+ws.nAS+1:end]
    MMAI = prob.MM[ws.AS,ws.IS]
    μTH .= @view prob.d[1:end-1,ws.IS]; mul!(μTH,λTH,MMAI,1,-1)
    μC .=  @view prob.d[end,ws.IS]; mul!(μC,MMAI',λC,1,1)
    return false, false
end

function compute_λ_and_μ(ws,prob::MPQP,opts)
    ws.nAS < prob.rank_defficiency && return false,true  # Trivially degenerate

    Q,R = qr(prob.A[ws.AS,:]')
    if(size(R,1) < ws.nAS || any(abs(R[i,i]) <1e-12 for i in 1:ws.nAS))
        return true, false # LICQ broken -> explore up, do not explore down
    end
    R = UpperTriangular(R)
    Q *= 1.0 # (to get explicit representation of Q) TODO: do inplace with lmul!
    Y,Z = Q[:,1:ws.nAS], Q[:,ws.nAS+1:end]

    Hr = Z'*prob.H*Z
    C = cholesky(Hr,check=false)
    if(!issuccess(C))
        return false,true# do not explore up, explore down
    end
    Rh = UpperTriangular(C.factors)

    xy = prob.B[:,ws.AS]/R
    xz = -(xy*Y'*prob.H+prob.F)*Z; rdiv!(rdiv!(xz, Rh), adjoint(Rh))
    x = xy*Y'+xz*Z'
    rhs = (x*prob.H+prob.F)*Y

    λTH = @view ws.Ath[:,ws.m0+1:ws.m0+ws.nAS]
    λC = @view ws.bth[ws.m0+1:ws.m0+ws.nAS]
    λTH .= @view rhs[1:end-1,:]; rdiv!(λTH, adjoint(R))
    λC .= -@view rhs[end,:]; ldiv!(R, λC)

    μTH = @view ws.Ath[:,ws.m0+ws.nAS+1:end]
    μC = @view ws.bth[ws.m0+ws.nAS+1:end]
    AIt = prob.A[ws.IS,:]
    μTH .= @view prob.B[1:end-1,ws.IS]; mul!(μTH,view(x,1:prob.n_theta,:),AIt',1,-1)
    μC .=  @view prob.B[end,ws.IS]; mul!(μC,AIt,x[end,:],-1,1)

    ws.x = x # Save x for later

    return false,false
end
## Explore subsets 
function explore_subsets(as,ws,id_cands,S)
    UIntX = typeof(as)
    for i in id_cands
        mask = UIntX(1)<<(i-1)
        as&mask ==  0 && continue
        as_new = as&~mask
        if as_new ∉ ws.explored
            push!(ws.explored,as_new) # Mark as explored
            push!(S,as_new)
        end
    end
end
## Explore supersets 
function explore_supersets(as,ws,id_cands,S,bounds_table)
    UIntX = typeof(as)
    for i in id_cands
        mask = UIntX(1)<<(i-1);
        as&(mask|(1<<(bounds_table[i]-1))) != 0 && continue
        #as&mask != 0 && continue
        as_new = as|mask
        if as_new ∉ ws.explored
            push!(ws.explored,as_new) # Mark as explored
            push!(S,as_new) # Put on stack
        end
    end
end

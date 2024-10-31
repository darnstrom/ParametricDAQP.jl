## mpsolve 
# Compute the explicit solution to multi-parameteric QP
function mpsolve(mpQP,Θ;opts=Settings(), AS0 = nothing) # bounds_table as option
    mpLDP = MPLDP(mpQP)
    mpLDP, Θ, tf = normalize_parameters(mpLDP,Θ)
    if(isnothing(AS0))
        AS0 = compute_AS0(mpLDP,Θ)
    end
    F,info =  mpdaqp_explicit(mpLDP,Θ,AS0;opts)
    return F, merge(info,(scaling=tf.scaling,translation=tf.center))
end
## Method based on combinatorial adjacency 
function mpdaqp_explicit(prob,Θ,AS0;opts = Settings())
    time_limit = opts.time_limit*1e9;
    status = :Solved
    t0  = time_ns()
    nth,m = size(prob.d);

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
            explore_supersets(pop!(ws.Sdown),ws,prob,m,ws.S)
        end
        while(length(ws.S) < opts.chunk_size && !isempty(ws.Sup)) # ... then try to move up
            explore_subsets(pop!(ws.Sup),ws,prob,m,ws.S)
        end

        opts.verbose>0 && print_ws(ws,(j+=1))

        # Process pending AS
        while(!isempty(ws.S))
            as = pop!(ws.S);
            region,LICQ_broken,down= isoptimal(as,ws,prob,opts)
            if(!isnothing(region))
                push!(ws.F,region)
                down && push!(ws.Sdown,as)
                push!(ws.Sup,as)
            end
            LICQ_broken && push!(ws.Sup,as)
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

    if(opts.factorization == :qr)
        R = (ws.nAS > 0) ? UpperTriangular(qr(prob.M[ws.AS,:]').R) : UpperTriangular(zeros(0,0))
        if(any(abs(R[i,i]) <1e-12 for i in 1:ws.nAS))
            return nothing,true,false # LICQ broken
        end
    else
        C = cholesky(prob.MM[ws.AS,ws.AS],check=false)
        if(!issuccess(C))
            return nothing,true,false # LICQ broken
        end
        R = UpperTriangular(C.factors)
    end

    # Compute λ
    λTH = @view ws.Ath[:,ws.m0+1:ws.m0+ws.nAS]
    λC = @view ws.bth[ws.m0+1:ws.m0+ws.nAS]
    compute_λ(prob,ws.AS,ws.IS,R,λTH,λC)

    # Compute μ
    μTH = @view ws.Ath[:,ws.m0+ws.nAS+1:end]
    μC = @view ws.bth[ws.m0+ws.nAS+1:end]
    compute_μ(prob,ws.AS,ws.IS,λTH,λC,μTH,μC)

    # Reset feasibility workspace
    reset_workspace(ws) 
    # Solve primal+dual feasibility problem 
    normalize_model(ws;eps_gap=opts.eps_gap,eps_zero=opts.eps_zero) || return nothing,false,false;

    ws.nLPs+=1
    if isfeasible(ws.DAQP_workspace, ws.m, 0)
        return extract_CR(ws,prob,λTH,λC)...,(ws.nAS ≤ prob.n)
    else
        return nothing,false,false
    end
end
## Update optimization model
function normalize_model(ws;eps_gap=1e-6,eps_zero=1e-12)
    norm_factor=0.0;
    # add λ ≥ 0 to feasibility model
    for i in (1+ws.m0: ws.m0+length(ws.AS))
        norm_factor = norm(view(ws.Ath,:,i),2);
        ws.norm_factors[i-ws.m0] = norm_factor
        if(norm_factor>eps_zero)
            rdiv!(view(ws.Ath,:,i),norm_factor)
            ws.bth[i]/=norm_factor
            ws.bth[i]-=eps_gap
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
function compute_λ(prob,AS,IS,R,λTH,λC)
    λTH .= @view prob.d[1:end-1,AS]; rdiv!(rdiv!(λTH, R), adjoint(R))
    λC .= -@view prob.d[end,AS]; ldiv!(R, ldiv!(adjoint(R), λC))
end
function compute_μ(prob,AS,IS,λTH,λC,μTH,μC)
    MMAI = prob.MM[AS,IS]
    μTH .= @view prob.d[1:end-1,IS]; mul!(μTH,λTH,MMAI,1,-1)
    μC .=  @view prob.d[end,IS]; mul!(μC,MMAI',λC,1,1)
end
## Explore subsets 
function explore_subsets(as,ws,prob,m,S)
    UIntX = typeof(as)
    for i in 1:m 
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
function explore_supersets(as,ws,prob,m,S)
    UIntX = typeof(as)
    for i in 1:m 
        mask = UIntX(1)<<(i-1);
        as&(mask|(1<<(prob.bounds_table[i]-1))) != 0 && continue
        #as&mask != 0 && continue
        as_new = as|mask
        if as_new ∉ ws.explored
            push!(ws.explored,as_new) # Mark as explored
            push!(S,as_new) # Put on stack
        end
    end
end

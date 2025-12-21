## Setup MPLDP, MPQP or MPVI
function setup_mpp(mpp;normalize=true, fix_ids=Int[],fix_vals=zeros(0))
    # If already setup, just return
    typeof(mpp) <: Union{MPLDP,MPQP,MPVI} && return deepcopy(mpp)

    if hasproperty(mpp,:f_theta)
        f_theta = mpp.f_theta 
    elseif hasproperty(mpp,:F)
        f_theta = mpp.F
    end

    if hasproperty(mpp,:W)
        W = mpp.W
    elseif hasproperty(mpp,:B)
        W = mpp.B
    end

    eq_ids = hasproperty(mpp,:eq_ids) && !isnothing(mpp.eq_ids) ? mpp.eq_ids : Int[]
    if hasproperty(mpp,:sense)
        eq_ids = eq_ids ∪ findall(mpp.sense.&DAQPBase.EQUALITY.!=0)
    elseif hasproperty(mpp,:senses)
        eq_ids = eq_ids ∪ findall(mpp.senses.&DAQPBase.EQUALITY.!=0)
    end

    n, nth = size(f_theta) 
    m = length(mpp.b)

    bnd_tbl = (haskey(mpp,:bounds_table) && !isnothing(mpp.bounds_table)) ? mpp.bounds_table : collect(1:m)
    out_inds = haskey(mpp,:out_inds) && !isnothing(mpp.out_inds) ? mpp.out_inds : collect(1:n)

    zlims =  get_lims(mpp.A,mpp.b,W,out_inds)


    free_ids = setdiff(1:nth,fix_ids)
    nth = length(free_ids) 

    F = [f_theta[:,free_ids] mpp.f+f_theta[:,fix_ids]*fix_vals]
    B = [W[:,free_ids] mpp.b+W[:,fix_ids]*fix_vals]

    # If H not symmetric => MPVI
    if(!isapprox(mpp.H, mpp.H', rtol=1e-9))
        AHinv = mpp.A / mpp.H
        HinvAt = mpp.A / mpp.H'
        AHinvA = HinvAt*mpp.A'
        #AHinvA = mpp.A*AHinv'
        HinvF = Matrix(-(mpp.H\F)')
        d = Matrix((B-mpp.A*HinvF')')
        return MPVI(mpp.H, Matrix(F'),mpp.A, Matrix(B'),HinvAt,AHinvA,HinvF,d,
                    nth,n,bnd_tbl,ones(m),eq_ids,out_inds,zlims)
    end

    # H is symmetric -> Either MPQP or MPLDP
    H = (mpp.H+mpp.H')/2
    R = cholesky(H, check=false)
    if(!issuccess(R)) # Cannot formulate as an LDP => return MPQP
        return MPQP(H,Matrix(F'),mpp.A,Matrix(B'),nth,n,bnd_tbl,ones(m),eq_ids,n-rank(H),out_inds,zlims)
    end

    # H is symmetric + positive definite => return MPLDP
    M = mpp.A/R.U
    V = (R.L)\F
    d = Matrix((B + M*V)')# Col. major...

    norm_factors = ones(m)
    if(normalize)
        norm_factor = 0
        for i in 1:m
            norm_factor = norm(view(M,i,:),2) 
            if(norm_factor > 0)
                M[i,:]./=norm_factor
                d[:,i]./=norm_factor
            end
            norm_factors[i]= norm_factor
        end
    end

    MRt = M/R.L
    RinvV = -(R.U\V)'
    MRt = MRt[:,out_inds]
    RinvV = RinvV[:,out_inds]

    return MPLDP(M*M', M, MRt, RinvV, d, nth, n, bnd_tbl, norm_factors, eq_ids,zlims)
end

## Normalize parameters to -1 ≤ θ ≤ 1
function normalize_parameters(prob::Union{MPLDP,MPVI},Θ;eps_zero=1e-12)
    if(isempty(Θ.ub)) # assume already normalized
        return prob,Θ,(center=0,scaling = 1)
    end
    nth = length(Θ.lb);
    center = (Θ.lb+Θ.ub)/2;
    norm_factors = (Θ.ub-Θ.lb)/2;
    prob.d[end:end,:] += center'*prob.d[1:end-1,:];
    prob.HinvF[end:end,:] += center'*prob.HinvF[1:end-1,:];
    for i in 1:nth
        prob.d[i,:] *= norm_factors[i];
        prob.HinvF[i,:] *= norm_factors[i];
    end

    Ath = haskey(Θ,:A) ? Θ.A : zeros(nth,0);
    bth = haskey(Θ,:b) ? Θ.b : zeros(0);
    A,b = normalize_region(Ath,bth,norm_factors,center)
    Θ =(A=A, b=b, lb=-ones(nth), ub=ones(nth));
    return prob, Θ,(center=center, scaling = 1 ./ norm_factors)
end
function normalize_parameters(prob::MPQP,Θ)
    if(isempty(Θ.ub)) # assume already normalized
        return prob,Θ,(center=0,scaling = 1)
    end
    nth = length(Θ.lb);
    center = (Θ.lb+Θ.ub)/2;
    norm_factors = (Θ.ub-Θ.lb)/2;
    prob.B[end:end,:] += center'*prob.B[1:end-1,:];
    prob.F[end:end,:] += center'*prob.F[1:end-1,:];
    for i in 1:nth
        prob.B[i,:] *= norm_factors[i];
        prob.F[i,:] *= norm_factors[i];
    end

    Ath = haskey(Θ,:A) ? Θ.A : zeros(nth,0);
    bth = haskey(Θ,:b) ? Θ.b : zeros(0);
    A,b = normalize_region(Ath,bth,norm_factors,center)
    Θ =(A=A,b=b, lb = -ones(nth), ub = ones(nth));
    return prob, Θ,(center=center, scaling = 1 ./ norm_factors)
end
function normalize_region(Ath,bth,norm_factors,center; eps_zero=1e-12)
    # Normalize A to box -1 ≤ θ ≤ 1
    A = norm_factors.*Ath
    b = bth-(center'*Ath)[:]

    # Normalize A (make half planes have norm = 1)
    is_nonzero = falses(length(b))
    for i in 1:length(b)
        norm_factor = norm(view(A,:,i),2);
        if(norm_factor>eps_zero)
            rdiv!(view(A,:,i),norm_factor)
            b[i]/=norm_factor
            is_nonzero[i] = true
        else
            (b[i]<-eps_zero) && return nothing,nothing,nothing # trivially infeasible
        end
    end
    return A[:,is_nonzero], b[is_nonzero]
end

## Denormalize parameters
function denormalize(v::AbstractVector,scaling,translation)
    return Diagonal(1 ./ scaling)*v+translation
end
function denormalize(A,b,scaling,translation)
    An = Diagonal(scaling)*A 
    bn = b+An'*translation
    return An,bn 
end
function denormalize(F::AbstractMatrix,scaling,translation;hps=false)
    Fn = Diagonal([scaling;1])*F 
    if hps
        Fn[end,:] += (translation'*Fn[1:end-1,:])[:]
    else
        Fn[end,:] -= (translation'*Fn[1:end-1,:])[:]
    end
    return Fn
end
function denormalize(cr::CriticalRegion,scaling,translation)
    if !isempty(cr.Ath)
        An,bn = denormalize(cr.Ath,cr.bth,scaling,translation)
    else 
        An,bn = zeros(0,0),zeros(0)
    end
    xn = !isempty(cr.z) ? denormalize(cr.z,scaling,translation) : cr.z 
    lamn = !isempty(cr.lam) ? denormalize(cr.lam,scaling,translation) : cr.lam
    thn = !isempty(cr.th) ? denormalize(cr.th,scaling,translation) : zeros(0)
    return CriticalRegion(cr.AS,An,bn,xn,lamn,thn)
end


## Reset workspace
function reset_workspace(ws) 
    ws.sense[1:ws.m].=0
    ws.m = ws.m0 
    p = ws.DAQP_workspace
    ccall((:reset_daqp_workspace,DAQPBase.libdaqp),Cvoid,(Ptr{Cvoid},),p);
end

## Setup workspace
# If ub/lb is not set, assume -1 ≤ θ ≤ 1 
function setup_workspace(Θ,n_constr;opts=Settings())::Workspace
    nth,n_general = size(Θ.A);
    m0 = isempty(Θ.ub) ? n_general : 2*nth+n_general;
    m_max = m0+n_constr
    A = Array{Float64}(undef,nth,m_max);
    b = Array{Float64}(undef,m_max);
    blower = fill(-1e30,m_max);
    # Setup Ath θ ≤ bth 
    if(!isempty(Θ.ub))
        A[:,1:nth] .= I(nth); 
        A[:,nth+1:2*nth] .= -I(nth); 
        b[1:nth] .= Θ.ub 
        b[nth+1:2*nth] .= -Θ.lb 
        A[:,2*nth+1:2*nth+n_general].=Θ.A;
        b[2*nth+1:2*nth+n_general].=Θ.b;
    else
        A[:,1:n_general].=Θ.A;
        b[1:n_general].=Θ.b;
    end

    UIntX = UInt
    if(n_constr <= 32)
        UIntX = UInt32
    elseif(n_constr <= 64)
        UIntX = UInt64
    elseif(n_constr <= 128)
        UIntX = UInt128
    elseif(n_constr <= 256)
        UIntX = UInt256
    elseif(n_constr <= 512)
        UIntX = UInt512
    elseif(n_constr <= 1024)
        UIntX = UInt1024
    else
        @error("Currently, 1024 is the maximum number of supported constraints")
    end
    # TODO: possible to use external package for larger UInt... 
    # 128 constraints is, however, probably above what is tractable anyways 

    #Create C workspace
    p=DAQPBase.setup_c_workspace(nth);
    # Set fval_bound to maximal radius for early termination
    # (the region is contained in a ball with this radius)
    max_radius =  isempty(Θ.ub) ? nth : nth*(maximum(Θ.ub)^2)/2;
    ws = Workspace{UIntX}(A,b,blower,zeros(Cint,m_max),0,m0,p,falses(0,0),0, 
                       UIntX[], UIntX[], UIntX[], CriticalRegion[],Set{UIntX}(),opts,
                       falses(n_constr),falses(n_constr),0, zeros(n_constr),zeros(0,0));
    DAQPBase.init_c_workspace_ldp(p,ws.Ath,ws.bth,ws.bth_lower,ws.sense;max_radius)
    settings(ws.DAQP_workspace,opts.daqp_settings)
    return ws 
end
## convert AS to unsigned int 
function as2uint(AS::Vector{Int},UIntX)
    as = UIntX(0)
    for i in AS 
        as |= UIntX(1) << (i-1)
    end
    return as
end

## extract AS and IS from unsigned int
function uint2as(u, ws, bounds_table)
    ws.AS .= false
    ws.IS .= true
    for i = 1:length(ws.AS)
        (u == 0) && break # the rest is inactive
        if(u&1==1)
            ws.AS[i] = true
            ws.IS[i] = false
            #ws.IS[i] = false;ws.IS[bounds_table[i]]=false;
        end
        u >>= 1
    end
    ws.nAS = sum(ws.AS)
end

## Extract Critical region 
function extract_CR(ws,prob)
    if(ws.opts.postcheck_rank) 
        islowrank(prob,ws) && return nothing
    end
    if(ws.opts.store_points)
        dws = unsafe_load(Ptr{DAQPBase.Workspace}(ws.DAQP_workspace))
        θ = copy(unsafe_wrap(Vector{Float64}, dws.u, dws.n, own=false))
    else
        θ = zeros(0)
    end
    if(ws.opts.lowdim_tol > 0)
        rhs_offset = ws.opts.lowdim_tol + 1e-6 # + 1e-6 to account for tolerance in DAQP
        ccall((:deactivate_constraints,DAQPBase.libdaqp),Cvoid,(Ptr{Cvoid},),ws.DAQP_workspace);
        ccall((:reset_daqp_workspace,DAQPBase.libdaqp),Cvoid,(Ptr{Cvoid},),ws.DAQP_workspace);
        ws.bth[1:ws.m].-=rhs_offset; # Shrink region 
        is_feasible = isfeasible(ws.DAQP_workspace, ws.m, 0)
        ws.bth[1:ws.m].+=rhs_offset; # Restore rhs
        !is_feasible &&  return nothing
    end

    AS = ws.opts.store_AS || ws.opts.store_regions ? findall(ws.AS) : Int64[]
    # Extract regions/solution
    if(ws.opts.store_regions)
        if(ws.opts.remove_redundant)
            ccall((:deactivate_constraints,DAQPBase.libdaqp),Cvoid,(Ptr{Cvoid},),ws.DAQP_workspace);
            Ath,bth = minrep(ws.DAQP_workspace); 
        else
            Ath,bth = ws.Ath[:,1:ws.m],ws.bth[1:ws.m];
        end
        x,λ = extract_solution(AS,prob,ws)
    else
        Ath,bth = zeros(prob.n_theta,0), zeros(0)
        x,λ = zeros(0,0),zeros(0,0);
    end

    return CriticalRegion(AS,Ath,bth,x,λ,θ)
end
## Check rank
islowrank(prob::MPLDP,ws) = rank(view(prob.AHinvA,ws.AS,ws.AS))<ws.nAS
islowrank(prob::MPQP,ws) = false # TODO
islowrank(prob::MPVI, ws) = false
## Extract solution
function extract_solution(AS,prob::Union{MPLDP,MPVI},ws)
    λ = [ws.Ath[:,ws.m0+1:ws.m0+ws.nAS];-ws.bth[ws.m0+1:ws.m0+ws.nAS]']
    # Renormalize dual variable from half-plane normalization
    for i in 1:ws.nAS
        rmul!(view(λ,:,i),ws.norm_factors[i])
    end

    # Compute primal solution x[1:end-1,:]' θ + x[end,:]
    x = copy(prob.HinvF); mul!(x,λ,prob.HinvAt[ws.AS,:],1,1)

    # Renormalize dual variable from LDP transform
    if ws.opts.store_dual && ws.opts.store_AS
        for (i,ASi) in enumerate(AS)
            rdiv!(view(λ,:,i),-prob.norm_factors[ASi])
        end
    end
    λ = ws.opts.store_dual ? λ : zeros(0,0)
    return x,λ
end
function extract_solution(AS,prob::MPQP,ws)
    if ws.opts.store_dual
        λ = [ws.Ath[:,ws.m0+1:ws.m0+ws.nAS];-ws.bth[ws.m0+1:ws.m0+ws.nAS]']
    else
        λ = zeros(0,0)
    end
    x = ws.x[:,prob.out_inds]
    return x,λ
end
## Compute AS0 
function compute_AS0(mpLDP::MPLDP,Θ)
    # Center in box is zero -> dtot = d[end,:]
    senses = zeros(Cint,size(mpLDP.d,2)+length(Θ.b));
    senses[mpLDP.eq_ids] .= DAQPBase.EQUALITY

    if(isempty(Θ.A) || all(Θ.b .>= 0)) # Check if origin is feasible
        _,_,exitflag,info= DAQPBase.quadprog(zeros(0,0),zeros(0),mpLDP.M,mpLDP.d[end,:],Float64[], senses);
        exitflag == 1 && return mpLDP.eq_ids ∪ findall(abs.(info.λ).> 0)
    end

    # Solve lifted feasibility problem in (x,θ)-space to find initial point 
    Alift = [-mpLDP.d[1:end-1,:]' mpLDP.M; Θ.A' zeros(length(Θ.b),mpLDP.n)]
    x,_,exitflag,info= DAQPBase.quadprog(zeros(0,0),zeros(0),Alift,[mpLDP.d[end,:];Θ.b],Float64[],senses);
    if exitflag != 1
        @warn "There is no parameter that makes the problem feasible"
        return nothing
    end
    θ = x[1:mpLDP.n_theta]
    _,_,exitflag,info= DAQPBase.quadprog(zeros(0,0),zeros(0),mpLDP.M,mpLDP.d'*[θ;1],Float64[],senses);
    return mpLDP.eq_ids ∪ findall(abs.(info.λ).> 0)
end
function compute_AS0(mpQP::MPQP,Θ)
    # Center in box is zero -> dtot = d[end,:]
    # TODO: set eps_prox ≠ 0
    senses = zeros(Cint,size(mpQP.B,2)+length(Θ.b));
    senses[mpQP.eq_ids] .= DAQPBase.EQUALITY
    d = DAQPBase.Model();
    DAQPBase.settings(d,Dict(:eps_prox=>1e-6)) # Since the Hessian is singular
    if(isempty(Θ.A) || all(Θ.b .>= 0)) # if origin is feasible
        DAQPBase.setup(d,mpQP.H,mpQP.F[end,:],mpQP.A,mpQP.B[end,:],Float64[], senses);
        x,fval,exitflag,info = DAQPBase.solve(d);
        exitflag == 1 && return findall(abs.(info.λ).> 0)
    end

    Alift = [-mpQP.B[1:end-1,:]' mpQP.A; Θ.A' zeros(length(Θ.b),mpQP.n)]
    # Solve lifted feasibility problem in (x,θ)-space to find initial point 
    x,_,exitflag,info= DAQPBase.quadprog(zeros(0,0),zeros(0),Alift,[mpQP.B[end,:];Θ.b],Float64[],senses);
    if exitflag != 1
        @warn "There is no parameter that makes the problem feasible"
        return nothing
    end
    θ = x[1:mpQP.n_theta]
    d = DAQPBase.Model();
    DAQPBase.settings(d,Dict(:eps_prox=>1e-6)) # Since the Hessian is singular
    DAQPBase.setup(d,mpQP.H,mpQP.F'*[θ;1],mpQP.A,mpQP.B'*[θ;1],Float64[],senses);
    x,fval,exitflag,info = DAQPBase.solve(d);
    return findall(abs.(info.λ).> 0)
end
function compute_AS0(mpQP::MPVI, Θ)
    # Solve lifted feasibility problem in (x,θ)-space to find initial point 
    senses = zeros(Cint, size(mpQP.B, 2)) # Vector that sets the type of constraints. 0 indicates inequality
    nth = size(Θ.A, 1)
    nx = size(mpQP.A, 2)
    n_constr_Θ = 0
    :A in keys(Θ) && (n_constr_Θ += size(Θ.A, 2))
    :ub in keys(Θ) && (n_constr_Θ += nth)
    :lb in keys(Θ) && (n_constr_Θ += nth)
    A_Θ = [Θ.A'; Matrix{Float64}(I(nth)); -Matrix{Float64}(I(nth))]
    b_Θ = [Θ.b; Θ.ub; -Θ.lb]
    A_lifted = [-mpQP.B[1:end-1, :]' mpQP.A; A_Θ zeros(n_constr_Θ, nx)]
    B_lifted = [mpQP.B[end, :]; b_Θ]
    x, _, exitflag, info = DAQPBase.quadprog(zeros(0,0), zeros(0), A_lifted, B_lifted)
    if exitflag != 1
        @warn "[compute_AS0] There is no parameter in the parameter set that makes the problem feasible"
        return nothing
    end
    θ = x[1:mpQP.n_theta]
    # Solve problem for the given θ
    tol = 10^(-4)
    x, r, solved = AVIsolve(mpQP.H, mpQP.F' * [θ; 1], mpQP.A, mpQP.B' * [θ; 1], max_iter=10^6, stepsize=0.9, tol=10^(-5))
    (solved == :MaximumIterationsReached) && @warn "[compute_AS0] Using inaccurate AVI solution"
    active_set = findall(mpQP.A * x .>= mpQP.B' * [θ; 1] .- tol)
    return active_set
end
## Get CRs 
function get_critical_regions(sol::Solution)
    return [denormalize(cr,sol.scaling,sol.translation) for cr in sol.CRs]
end
## Get unsaturated regions
function get_unsaturated(CRs::Vector{CriticalRegion};tol = 1e-5)
    fbs, reg2fb = get_feedbacks(CRs)
    Nr = length(CRs)
    unsat_ids = Int[]
    for (k,f) in enumerate(fbs)
        if  norm(f[1:end-1,:],Inf) > tol
            append!(unsat_ids,findall(reg2fb[i] == k for i in 1:Nr))
        end
    end
    return CRs[unsat_ids]
end
## Get limits of out variables
function get_lims(A,b,W,out_inds)
    zlims = [-1e30*ones(length(out_inds)) 1e30*ones(length(out_inds))]
    for (k,id) in enumerate(out_inds)
        for (i,a) = enumerate(eachrow(A))
            if(a[id] == 1) &&  sum(!=(0),a) == 1 && iszero(W[i,:])
                zlims[k,2] = min(b[i],zlims[k,2])
            elseif(a[id] == -1) &&  sum(!=(0),a) == 1
                zlims[k,1] = max(-b[i],zlims[k,1])
            end
        end
    end
    return zlims
end

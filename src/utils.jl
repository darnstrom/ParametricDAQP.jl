## Setup MPLDP from MPQP
function MPLDP(mpQP;normalize=true)
    n, nth = size(mpQP.f_theta) 
    m = length(mpQP.b)

    R = cholesky((mpQP.H+mpQP.H')/2)
    M = mpQP.A/R.U
    V = (R.L)\[mpQP.f_theta mpQP.f]
    d = Matrix(([mpQP.W mpQP.b] + M*V)')# Col. major...

    if(normalize)
        norm_factor = 0
        for i in 1:m
            norm_factor = norm(view(M,i,:),2) 
            M[i,:]./=norm_factor
            d[:,i]./=norm_factor
        end
    end

    # unconstrained solution given by vTH θ + vC
    vTH = -(R.U\V[:,1:end-1])'
    vC  = -(R.U\V[:,end])

    bnd_tbl = haskey(mpQP,:bounds_table) ? mpQP.bounds_table : collect(1:m)
    haskey(mpQP,:out_inds) && (MRt = MRt[:,mpQP.out_inds])

    return MPLDP(M*M', M, (M/R.L), vTH, vC, d, nth, n, bnd_tbl)
end

## Normalize parameters to -1 ≤ θ ≤ 1
function normalize_parameters(prob,Θ)
    if(isempty(Θ.ub)) # assume already normalized
        return prob,Θ,(center=0,scaling = 1)
    end
    nth = length(Θ.lb);
    center = (Θ.lb+Θ.ub)/2;
    norm_factors = (Θ.ub-Θ.lb)/2;
    prob.d[end:end,:] += center'*prob.d[1:end-1,:];
    for i in 1:nth
        prob.d[i,:] *= norm_factors[i];
    end

    Ath = haskey(Θ,:A) ? Θ.A : zeros(nth,0);
    bth = haskey(Θ,:b) ? Θ.b : zeros(0);
    Θ =(A=norm_factors.*Ath, b = bth-(center'*Ath)[:], # TODO: verify
        lb = -ones(nth), ub = ones(nth));
    return prob, Θ,(center=center, scaling = 1 ./ norm_factors)
end

## Reset workspace
function reset_workspace(ws) 
    ws.sense[1:ws.m].=0
    ws.m = ws.m0 
    p = ws.DAQP_workspace
    ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);
end

## Setup workspace
# If ub/lb is not set, assume -1 ≤ θ ≤ 1 
function setup_workspace(Θ,n_constr;opts=EMPCSettings())::EMPCWorkspace
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
    else
        @error("Currently, 128 is the maximum number of supported constraints")
    end
    # TODO: possible to use external package for larger UInt... 
    # 128 constraints is, however, probably above what is tractable anyways 

    #Create C workspace
    p=DAQP.setup_c_workspace(nth);
    # Set fval_bound to maximal radius for early termination
    # (the region is contained in a ball with this radius)
    max_radius =  isempty(Θ.ub) ? nth : nth*(maximum(Θ.ub)^2); 
    ws = EMPCWorkspace{UIntX}(A,b,blower,zeros(Cint,m_max),0,m0,p,falses(0,0),0, 
                       UIntX[], UIntX[], UIntX[], CriticalRegion[],Set{UIntX}(),opts,
                       falses(n_constr),falses(n_constr),0, zeros(n_constr));
    DAQP.init_c_workspace_ldp(p,ws.Ath,ws.bth,ws.bth_lower,ws.sense;max_radius)
    return ws 
end
## convert AS to unsigned int 
function as2uint(AS::Vector{Int},UIntX)
    as = UIntX(0)
    for i in AS 
        as |= 1 << (i-1)
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
function extract_CR(ws,prob,λTH,λC)
    if(ws.opts.postcheck_rank && rank(view(prob.MM,ws.AS,ws.AS))<ws.nAS)
        return nothing,true
    end
    if(ws.opts.lowdim_tol > 0)
        ws.sense[1:ws.m].=0
        ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),ws.DAQP_workspace);
        ws.bth[1:ws.m].-=ws.opts.lowdim_tol; # Shrink region 
        if !isfeasible(ws.DAQP_workspace, ws.m, 0) # Check if region is narrow and, hence, should be removed
            return nothing,true
        end
    end

    AS = ws.opts.store_AS ? findall(ws.AS) : Int64[]

    if(ws.opts.store_points)
        dws = unsafe_load(Ptr{DAQP.Workspace}(ws.DAQP_workspace))
        θ = copy(unsafe_wrap(Vector{Float64}, dws.u, dws.n, own=false))
    else
        θ = zeros(0)
    end
    # Extract regions/solution
    if(ws.opts.store_regions)
        if(ws.opts.remove_redundant)
            Ath,bth = minrep(ws.DAQP_workspace); 
        else
            Ath,bth = ws.Ath[:,1:ws.m],ws.bth[1:ws.m];
        end
        # Compute primal solution xTH θ + xC 
        MRt = ws.norm_factors[1:ws.nAS].*prob.MRt[ws.AS,:]
        xTH = copy(prob.vTH); mul!(xTH,λTH,MRt,1,1) 
        λC .+= ws.opts.eps_gap; # Shift back (ok since ws.Ath already copied 
        xC = copy(prob.vC); mul!(xC,MRt',λC,-1,1)
    else
        Ath,bth = zeros(prob.n_theta,0), zeros(0)
        xTH = zeros(0,0); xC = zeros(0)
    end

    return CriticalRegion(AS,Ath,bth,xTH,xC,θ),false
end
## Compute AS0 
function compute_AS0(mpQP,θ)
    f = mpQP.f[:,1]+mpQP.f_theta*θ
    bu = mpQP.b[:,1]+mpQP.W*θ
    bl = -1e30*ones(length(mpQP.b))
    senses= haskey(mpQP,:senses) ? mpQP.senses : zeros(Cint,length(mpQP.b));

    _,_,_,info= DAQP.quadprog(mpQP.H,f,mpQP.A ,bu ,bl, senses);
    return findall(abs.(info.λ).> 0)
end

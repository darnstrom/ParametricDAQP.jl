## Setup MPLDP from MPQP
function MPLDP(mpQP;normalize=true)
    if hasproperty(mpQP,:f_theta)
        f_theta = mpQP.f_theta 
    elseif hasproperty(mpQP,:F)
        f_theta = mpQP.F
    end

    if hasproperty(mpQP,:W)
        W = mpQP.W
    elseif hasproperty(mpQP,:B)
        W = mpQP.B
    end

    n, nth = size(f_theta) 
    m = length(mpQP.b)

    R = cholesky((mpQP.H+mpQP.H')/2)
    M = mpQP.A/R.U
    V = (R.L)\[f_theta mpQP.f]
    d = Matrix(([W mpQP.b] + M*V)')# Col. major...

    norm_factors = ones(m)
    if(normalize)
        norm_factor = 0
        for i in 1:m
            norm_factor = norm(view(M,i,:),2) 
            M[i,:]./=norm_factor
            d[:,i]./=norm_factor
            norm_factors[i]= norm_factor
        end
    end

    bnd_tbl = (haskey(mpQP,:bounds_table) && !isnothing(mpQP.bounds_table)) ? mpQP.bounds_table : collect(1:m)

    MRt = M/R.L
    RinvV = -(R.U\V)'
    if haskey(mpQP,:out_inds) && !isnothing(mpQP.out_inds)
        MRt = MRt[:,mpQP.out_inds]
        RinvV = RinvV[:,mpQP.out_inds]
    end

    return MPLDP(M*M', M, MRt, RinvV, d, nth, n, bnd_tbl, norm_factors)
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
        prob.RinvV[i,:] *= norm_factors[i];
    end

    Ath = haskey(Θ,:A) ? Θ.A : zeros(nth,0);
    bth = haskey(Θ,:b) ? Θ.b : zeros(0);
    Θ =(A=norm_factors.*Ath, b = bth-(center'*Ath)[:], # TODO: verify
        lb = -ones(nth), ub = ones(nth));
    return prob, Θ,(center=center, scaling = 1 ./ norm_factors)
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
    xn = !isempty(cr.x) ? denormalize(cr.x,scaling,translation) : zeros(0,0)
    lamn = !isempty(cr.lam) ? denormalize(cr.lam,scaling,translation) : zeros(0,0)
    thn = !isempty(cr.th) ? denormalize(cr.th,scaling,translation) : zeros(0)
    return CriticalRegion(cr.AS,An,bn,xn,lamn,thn)
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
    ws = Workspace{UIntX}(A,b,blower,zeros(Cint,m_max),0,m0,p,falses(0,0),0, 
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
    if(ws.opts.store_points)
        dws = unsafe_load(Ptr{DAQP.Workspace}(ws.DAQP_workspace))
        θ = copy(unsafe_wrap(Vector{Float64}, dws.u, dws.n, own=false))
    else
        θ = zeros(0)
    end
    if(ws.opts.lowdim_tol > 0)
        ws.sense[1:ws.m].=0
        ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),ws.DAQP_workspace);
        ws.bth[1:ws.m].-=ws.opts.lowdim_tol; # Shrink region 
        if !isfeasible(ws.DAQP_workspace, ws.m, 0) # Check if region is narrow and, hence, should be removed
            return nothing,true
        end
        ws.bth[1:ws.m].+=ws.opts.lowdim_tol; # Restore 
    end

    AS = ws.opts.store_AS ? findall(ws.AS) : Int64[]
    # Extract regions/solution
    if(ws.opts.store_regions)
        if(ws.opts.remove_redundant)
            Ath,bth = minrep(ws.DAQP_workspace); 
        else
            Ath,bth = ws.Ath[:,1:ws.m],ws.bth[1:ws.m];
        end
        # Renormalize dual variable from half-plane normalization
        λ = [λTH;-λC']
        for i in 1:ws.nAS
            rmul!(view(λ,:,i),ws.norm_factors[i])
        end

        # Compute primal solution x[1:end-1,:]' θ + x[end,:]
        x = copy(prob.RinvV); mul!(x,λ,prob.MRt[ws.AS,:],1,1)
 
        # Renormalize dual variable from LDP transform
        if ws.opts.store_dual & ws.opts.store_AS
            for (i,ASi) in enumerate(AS)
                rdiv!(view(λ,:,i),-prob.norm_factors[ASi])
            end
        end
    else
        Ath,bth = zeros(prob.n_theta,0), zeros(0)
        x,λ = zeros(0,0),zeros(0,0);
    end

    return CriticalRegion(AS,Ath,bth,x,λ,θ),false
end
## Compute AS0 
function compute_AS0(mpLDP,Θ)
    # Center in box is zero -> dtot = d[end,:]
    _,_,_,info= DAQP.quadprog(zeros(0,0),zeros(0),mpLDP.M,mpLDP.d[end,:]);
    return findall(abs.(info.λ).> 0)
    # TODO add backup if this fails
end
## Get CRs 
function get_critical_regions(sol::Solution)
    return [denormalize(cr,sol.scaling,sol.translation) for cr in sol.CRs]
end

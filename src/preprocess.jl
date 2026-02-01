function preprocess_eliminate_equalities(mpp)
    ncstr = length(mpp.bu);
    n_bounds = ncstr-size(mpp.A,1);
    # Expand A 
    A0 = [I(n_bounds) zeros(n_bounds,size(mpp.A,2)-n_bounds);mpp.A]

    # Get mask for equalities
    mask = mpp.senses .== DAQPBase.EQUALITY
    eq_ids = findall(mask) 
    ineq_ids = findall(.~mask)

    # Extract equalities
    Aeq = A0[eq_ids,:]
    beq = mpp.bu[eq_ids]
    Beq = mpp.B[eq_ids,:]
    senses_eq = mpp.senses[eq_ids]

    # Extract inequalities
    Aineq = A0[ineq_ids,:]
    bu,bl = mpp.bu[ineq_ids],mpp.bl[ineq_ids]
    Bineq = mpp.B[ineq_ids,:]
    senses_ineq = mpp.senses[ineq_ids]

    # Eliminate equalities 
    z0,zTH,Z = Aeq\beq, Aeq\Beq, nullspace(Aeq)
    f,F = Z'*(mpp.f+mpp.H*z0), Z'*(mpp.F+mpp.H*zTH) 
    H = Z'mpp.H*Z # TODO cholesky -> QR
    bl, bu =  bl - Aineq*z0, bu - Aineq*z0
    B, Aineq = Bineq-Aineq*zTH, Aineq*Z

    return (H=H,F=F,f=f,
            A=Aineq,B=B,bl=bl,bu=bu,senses=senses_ineq, 
            post_transform = (Z,z0,zTH))
end

function make_singlesided(mpp;single_soft=false, eliminate_equalities=false,
        explicit_soft= true, soft_weight=1e6)
    ncstr = length(mpp.bu);
    n_bounds = ncstr-size(mpp.A,1);
    A0 = [I(n_bounds) zeros(n_bounds,size(mpp.A,2)-n_bounds);mpp.A]

    if hasproperty(mpp,:senses)
        senses = mpp.senses
    elseif hasproperty(mpp,:sense)
        senses = mpp.sense
    else
        senses = zeros(Cint,ncstr)
    end

    if hasproperty(mpp,:W)
        B = [mpp.W;-mpp.W]
    elseif hasproperty(mpp,:B)
        B = [mpp.B;-mpp.B]
    end

    if hasproperty(mpp,:prio)
        prio = mpp.prio
    else
        prio = zeros(Cint,ncstr)
    end

    nth = size(B,2)

    mask = senses .== DAQPBase.EQUALITY
    eq_ids = findall(mask) 
    ineq_ids = findall(.~mask)
    me = sum(mask)
    mi = ncstr-me

    # Extract equalities
    Aeq = A0[eq_ids,:]
    beq = mpp.bu[eq_ids]
    Beq = mpp.B[eq_ids,:]
    senses_eq = mpp.senses[eq_ids]
    prio_eq = prio[eq_ids]

    # Extract inequalities
    Aineq = A0[ineq_ids,:]
    bu,bl = mpp.bu[ineq_ids],mpp.bl[ineq_ids]
    Bineq = mpp.B[ineq_ids,:]
    senses_ineq = mpp.senses[ineq_ids]
    prio_ineq = prio[ineq_ids]


    A = [Aeq;Aineq;-Aineq]
    b = [beq;bu;-bl]
    B = [Beq;Bineq;-Bineq]
    senses = [senses_eq;senses_ineq;senses_ineq]
    prio = [prio_eq;prio_ineq;prio_ineq]
    bounds_table = [1:me; me+mi+1:me+2*mi;me+1:me+mi]

    # === Objective ====
    H,f,= mpp.H, mpp.f

    H_theta = hasproperty(mpp,:H_theta) ? mpp.H_theta : zeros(nth,nth)

    if hasproperty(mpp,:f_theta)
        F = mpp.f_theta
    elseif hasproperty(mpp,:F)
        F = mpp.F
    end

    # Introduce slack variables 
    if(explicit_soft)
        soft_mask = ((senses[1:me+mi] .& DAQPBase.SOFT) .== DAQPBase.SOFT)
        if(any(soft_mask))
            soft_ids = findall(soft_mask)
            R = cholesky((mpp.H+mpp.H')/2)
            Ms= A0[soft_mask,:]/R.U
            norm_factors = [norm(view(Ms,i,:),2) for i in 1:size(Ms,1)]
            if(single_soft)
                nsoft = 1
                A = [A zeros(me+2mi,1)]
                A[soft_ids,end] .= -norm_factors
                A[bounds_table[soft_ids],end] .= -norm_factors
            else
                nsoft, n = length(soft_ids), size(A,2)
                A = [A zeros(me+2mi,nsoft)]
                A[soft_ids,n+1:end] = -diagm(norm_factors)
                A[bounds_table[soft_ids],n+1:end] = -diagm(norm_factors)
            end
            H = cat(H,soft_weight*I(nsoft),dims=(1,2))
            f = [f;zeros(nsoft)] 
            F = [F;zeros(nsoft,size(F,2))]
        end
    end

    # Prune possible Inf bounds
    rm_ids = findall(b[:] .>= 1e20)
    if(!isempty(rm_ids))
        bounds_table[bounds_table[rm_ids]] = bounds_table[rm_ids] # Make other bound point to itself
        # Correct bounds table 
        rm_offset, keep_ids = 1, Int[]
        for i in 1:(me+2*mi)
            if(rm_offset <=length(rm_ids) && i==rm_ids[rm_offset])
                rm_offset+=1
            else
                push!(keep_ids,i)
            end
        end
        A,b,B = A[keep_ids,:],b[keep_ids],B[keep_ids,:]
        senses,prio,bounds_table = senses[keep_ids],prio[keep_ids],bounds_table[keep_ids]
    end

    # TODO fix outputs

    return (H=H,f=f, H_theta = H_theta, F=F,
            A=Matrix{Float64}(A), b=b, B=B, senses=senses,
            bounds_table=bounds_table, prio=prio)
end

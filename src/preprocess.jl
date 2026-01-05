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
    me = sum(mask_eq)
    mi = ncstr-me

    Aeq = A0[eq_ids,:]
    beq = mpp.bu[eq_ids]
    Beq = B[eq_ids,:]
    senses_eq = senses[eq_ids]
    prio_eq = prio[eq_ids]

    Aineq = mpp.A[ineq_ids,:]
    bu,bl = mpp.bu[ineq_ids],mpp.bl[ineq_ids]
    Bineq = B[ineq_ids,:]
    senses_ineq = senses[ineq_ids]
    prio_ineq = prio[ineq_ids]

    A = [Aeq;Aineq;-Aineq]
    b = [beq;bu;-bl]
    B = [Beq;Bineq;-Bineq]
    senses = [senses_eq;senses_ineq;sense_ineq]
    prio = [prio_eq;prio_eq;prio_eq]
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
        for i in 1:2*ncstr
            if(i==rm_ids[rm_offset])
                rm_offset+=1
            else
                bounds_table[i] -= (rm_offset-1)
                push!(keep_ids,i)
            end
        end
        A,b,B = A[keep_ids,:],b[keep_ids],B[keep_ids,:]
        senses,prio,bounds_table = senses[keep_ids],prio[keep_ids],bounds_table[keep_ids]
    end

    # Eliminate equalities
    if(eliminate_equalities && me > 0)
        # TODO
    else
        out_transform = I
    end

    return (H=H,f=f, H_theta = H_theta, F=F,
            A=Matrix{Float64}(A), b=b, B=B, senses=senses,
            bounds_table=bounds_table, prio=prio)
end

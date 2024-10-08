function form_tree(F;δ=0.01,depth=3)
    nth = size(F[1].Ath,1)
    bbs = []
    # Compute bounding boxes for critical regions 
    for f in F
        push!(bbs,PolyDAQP.bounding_box(f.Ath,f.bth))
    end

    hs = Vector{Float64}[]
    for i = 1:nth
        push!(hs,([-1.0,1.0]))
    end

    for bb in bbs
        ub,lb = bb
        for i = 1:nth
            for h in (ub[i],lb[i])
                if(all(abs.(hs[i].-h).> δ))
                    push!(hs[i],h)
                end
            end
        end
    end
    sort!.(hs)
    

    Hs = []
    for i = 1:nth
        for h in hs[i]
            push!(Hs,(i,h))
        end
    end

    Ip, In = [],[] 
    for h in Hs
        i,val = h
        Iph,Inh = [],[]
        for (reg_id,bb) in enumerate(bbs)
            ub,lb = bb
            if(ub[i] >= val )
                push!(Iph,reg_id)
            end
            if(lb[i] <= val)
                push!(Inh,reg_id)
            end
        end
        push!(Ip,Iph)
        push!(In,Inh)
    end


    U = [(collect(1:length(bbs)),Int[])]
    tree = [];
    while !isempty(U)
        regs, H_ids = pop!(U)
        add_cand,val = 0,Inf;
        Ipc,Inc = [],[];
        for (k,h) in enumerate(Hs)
            # TODO: Maybe minimize on the number of "free constraints" instead?
            Ipk = regs ∩ Ip[k]
            (length(Ipk) > val) && continue 
            Ink = regs ∩ In[k]
            (length(Ink) > val) && continue 
            add_cand,val = k, max(length(Ipk),length(Ink))
            Ipc,Inc = Ipk,Ink
        end
        #TODO: Get true positive/negative partition by solving LP 

        if length(H_ids) == depth  
            push!(tree,(Ipc,[H_ids;add_cand]))
            push!(tree,(Inc,[H_ids;-add_cand]))
        else
            push!(U,(Ipc, [H_ids;add_cand]))
            push!(U,(Inc, [H_ids;-add_cand]))
        end
    end

    return tree 
end

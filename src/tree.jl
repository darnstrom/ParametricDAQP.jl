struct BinarySearchTree
    halfplanes::Matrix{Float64}
    feedbacks::Vector{Matrix{Float64}}
    hp_list::Vector{Int}
    jump_list::Vector{Int}
    depth::Int
    duals::Vector{Matrix{Float64}}
    clipping::Matrix{Float64}
end

function get_halfplanes(CRs;tol=1e-5)
    nreg = length(CRs)
    nreg == 0 && return nothing
    nth = size(CRs[1].Ath,1)
    hps = Vector{Float64}[]
    reg2hp = [Dict{Int,Int}() for _ in 1:nreg]

    for (reg_id,cr) in enumerate(CRs)
        for (i,a) = enumerate(eachcol(cr.Ath))
            # Disregard box bounds
            cr.bth[i] == 1 && sum(!=(0),a) == 1 && continue
            nz_id = findfirst(!=(0),a)
            isnothing(nz_id) && continue
            asign = sign(a[nz_id])
            hcand = asign*[a;cr.bth[i]]
            # Check if hcand already exists in hps
            new_hp = true
            for (j,h) in  enumerate(hps)
                if(all(isapprox(h[i],hcand[i],atol=tol,rtol=tol) for i in 1:nth+1))
                    reg2hp[reg_id][j] = asign
                    new_hp = false
                    break;
                end
            end
            if new_hp
                push!(hps, hcand)
                reg2hp[reg_id][length(hps)]=asign
            end
        end
    end
    return hcat(hps...),reg2hp
end

function get_feedbacks(CRs; tol=1e-5)
    isempty(CRs) && return nothing
    nth = size(CRs[1].Ath,1)
    Z,ids = Matrix{Float64}[], Int[]
    for cr in CRs 
        id = 0 
        for (i,z) in enumerate(Z)
            if all(isapprox.(cr.z,z,atol=tol, rtol=tol))
                id = i
                push!(ids,id)
                break
            end
        end
        if id == 0 # No duplicate
            push!(Z,cr.z)
            push!(ids,length(Z))
        end
    end
    return Z,ids
end

# TODO: Can be cut in half by using points in CR 
function classify_regions(CRs,hps, reg2hp, ws; reg_ids = nothing, hp_ids = nothing, branches = nothing)
    eps_gap=1e-6+1e-12
    reg_ids = isnothing(reg_ids) ?  (1:length(CRs)) : findall(reg_ids)
    isnothing(hp_ids) && (hp_ids = 1:size(hps,2))
    nth, nh = size(hps,1)-1, length(hp_ids)
    
    nR = length(CRs)
    nregs = [falses(nR) for _ in 1:nh]
    pregs = [falses(nR) for _ in 1:nh]


    nbr = isnothing(branches) ? 0 : length(branches)
    if !isnothing(branches)
        for i = 1:nbr
            (hid,hsign) = branches[i]
            @views ws.A[:,1+i] .= hsign.*hps[1:nth,hid]
            ws.b[1+i] = hsign*hps[end,hid]-eps_gap
        end
    end

    for i in reg_ids
        mi = length(CRs[i].bth)
        ws.A[:,1+nbr+1:1+nbr+mi] = CRs[i].Ath
        ws.b[1+nbr+1:1+nbr+mi] = CRs[i].bth .-eps_gap

        for (j,hj) in enumerate(hp_ids)
            # First check if the hp is a facet of the region
            asign = get(reg2hp[i], hj, nothing)
            if !isnothing(asign) # hp is a facet of the region
                if asign == 1
                    pregs[j][i] = true
                else
                    nregs[j][i] = true
                end
                continue
            end

            @views slack = isnothing(branches) ? hps[1:end-1,hj]'*CRs[i].th-hps[end,hj] : NaN

            # Negative
            @views ws.A[:,1] .= .-hps[1:nth,hj]
            @views ws.b[1] = -hps[end,hj]-eps_gap
            (slack > eps_gap || isfeasible(ws.p, 1+nbr+mi, 0)) && (nregs[j][i] = true)

            # Positive
            @views ws.A[:,1] = hps[1:nth,hj]
            @views ws.b[1] = hps[end,hj]-eps_gap
            (slack < -eps_gap || isfeasible(ws.p, 1+nbr+mi, 0)) && (pregs[j][i] = true)
        end
    end
    return nregs,pregs
end

function reduce_candidates(criteria,splits,ids)
    vals = [criteria(s) for s in splits]
    min_val = minimum(vals)
    min_ids = findall(==(min_val),vals)
    return min_val,splits[min_ids],ids[min_ids]
end

function get_split(CRs,hps,reg2hp,reg_ids,pregs,nregs,branches,criterions,ws,hp_ids_bit;balancing_level=0)
    fill!(hp_ids_bit,false)
    for i in findall(reg_ids)
        for hp in reg2hp[i]
            hp_ids_bit[first(hp)] = true
        end
    end
    for b in branches
        hp_ids_bit[first(b)] = false
    end
    hp_ids = findall(hp_ids_bit)

    #hp_ids = reduce(∪,Set(first.(reg2hp[i])) for i in findall(reg_ids));
    #hp_ids = collect(setdiff!(hp_ids,first(b) for b in branches))

    splits = [(reg_ids .* nregs[i], reg_ids .* pregs[i]) for i in hp_ids]
    isempty(splits) && return hp_ids, (falses(0),falses(0))

    # Use heuristic to find candidates
    if(length(branches) >= balancing_level || length(branches) == 0)
        min_val,splits,hp_ids = reduce_candidates(criterions[1],splits,hp_ids)
    else
        min_val = Inf
    end
    if length(branches) > 0 && min_val > 1# Compute the actual split
        splits = tuple.(classify_regions(CRs,hps,reg2hp,ws;reg_ids,hp_ids,branches)...)
        for c in criterions
            min_val,splits,hp_ids = reduce_candidates(c,splits,hp_ids)
            length(splits) == 1 && break
        end
    end
    hp_id = hp_ids[1]
    new_nregs, new_pregs = splits[1]
    return hp_ids[1], splits[1]
end

function setup_daqp_workspace(n_theta, nc_max = 250)
    p=DAQPBase.setup_c_workspace(n_theta);
    A = Array{Float64}(undef,n_theta,nc_max);
    b = Array{Float64}(undef,nc_max);
    bl =-1e30*ones(nc_max)
    sense = zeros(Cint,nc_max)

    DAQPBase.init_c_workspace_ldp(p,A,b,bl,sense;max_radius=n_theta/2)
    d_work = unsafe_load(Ptr{DAQPBase.Workspace}(p));
    unsafe_store!(Ptr{Cdouble}(d_work.settings+fieldoffset(DAQPBase.DAQPSettings,14)),1e-9); # sing_tol
    return (p=p,A=A,b=b,bl=bl,sense=sense)
end

function remove_redundant_hps(jump_list,hp_list,hps)
    hp_rows = findall(jump_list.!== 0)
    new_hp_ids = sort(unique(hp_list[hp_rows]))
    hp_old2new = Dict(id => k for (k,id) in enumerate(new_hp_ids))
    new_hps = hps[:,new_hp_ids]

    new_hp_list = copy(hp_list)
    for i in hp_rows
        new_hp_list[i] = hp_old2new[new_hp_list[i]]
    end
    return new_hps,new_hp_list
end

function get_duals(CRs,sol)
    fbs_dual = Matrix{Float64}[]
    for cr in CRs
        lam = zeros(size(cr.lam,1),length(sol.problem.norm_factors))
        lam[:,cr.AS] = cr.lam
        push!(fbs_dual,lam)
    end
    return [denormalize(f,sol.scaling,sol.translation) for f in fbs_dual]
end

function build_tree(sol::Solution; daqp_settings = nothing, verbose=1, max_reals=1e12,
        dual=false, bfs=true, clipping=false, balancing_level=0, hp_tol = 1e-5)
    if sol.status != :Solved
        verbose > 0 && @warn "Cannot build binary search tree. Solution status: $(sol.status)"
        return nothing
    end
    verbose > 0 && @info "Building binary search tree" 

    CRs = clipping ? get_unsaturated(sol.CRs) : sol.CRs

    # Get halfplanes and feedbacks
    hps,reg2hp = get_halfplanes(CRs;tol=hp_tol)
    fbs, fb_ids = get_feedbacks(CRs)

    nfbs,Nr = length(fbs),length(CRs)
    # Approximate number of real numbers
    n_reals = length(hps) + nfbs*length(fbs[1])
    dual && (n_reals += Nr*(sol.problem.n_theta+1)*length(sol.problem.norm_factors))
    if n_reals > max_reals
        verbose > 0 && @warn "Memory limit for real numbers (n_reals = $n_reals, max_reals = $max_reals) reached"
        return nothing
    end

    ws = setup_daqp_workspace(sol.problem.n_theta)

    # Do initial classification
    nregs,pregs = classify_regions(CRs,hps,reg2hp,ws)

    function get_extream_fbs(splits::Tuple{BitVector,BitVector}, fb_ids, seen_buffer, extreama)
        s_left, s_right = splits

        l_count = count_unique_ids(s_left, fb_ids, seen_buffer)
        r_count = count_unique_ids(s_right, fb_ids, seen_buffer)
        return extreama(l_count, r_count)
    end

    function count_unique_ids(bits, fb_ids, seen_buffer)
        fill!(seen_buffer, false)
        @inbounds for i in findall(bits)
            seen_buffer[fb_ids[i]] = true
        end
        return sum(seen_buffer)
    end

    seen_buffer = [false for _ in 1:length(fbs)]
    hp_ids_bit = falses(size(hps,2))

    main_criterion = (Nr == nfbs) ? s->max(sum.(s)...) : s->get_extream_fbs(s,fb_ids,seen_buffer,max)


    get_fbid = s->Set{Int}(fb_ids[s])
    criterions = dual ? [x->max(sum.(x)...)] :  [main_criterion,
                                                 s->sum(s[1] .& s[2]),
                                                 s->max(sum.(s)...),
                                                 s->get_extream_fbs(s,fb_ids,seen_buffer,min)]
    # Start exploration
    hp_list, jump_list = Int[0],Int[0]
    depth = 0
    tree_pop! = bfs ? popfirst! : pop!
    U = [(trues(length(CRs)),[],1)]
    while !isempty(U)
        reg_ids, branches, self_id = tree_pop!(U)
        depth = max(depth,length(branches))
        # Get halfplane to cut
        hp_id, (new_nregs, new_pregs) = get_split(CRs,hps,reg2hp,reg_ids,pregs,nregs,branches,criterions,ws,
                                                  hp_ids_bit; balancing_level)
        if isempty(hp_id) # Should never happen, but might due to numerics
            jump_list[self_id] = 0 # pointing at root node -> leaf
            hp_list[self_id] = first(get_fbid(reg_ids))
            @debug "Superfluous branch -> might be due to bad numerics"
            continue
        end
        #vals = [c((new_nregs,new_pregs)) for c in criterions]

        # Update tree for current node
        next_id = length(hp_list)+1
        hp_list[self_id] = hp_id
        jump_list[self_id] = next_id-self_id
        
        # Make room for the two new nodes that are spawned 
        push!(hp_list,0,0)
        push!(jump_list,0,0) 

        # Spawn new nodes
        for (new_regs,new_regs_comp,next, hp_sign) in [(new_nregs,new_pregs,next_id,-1), (new_pregs,new_nregs,next_id+1,1)]
            fb_cands = get_fbid(new_regs)
            if length(fb_cands) == 0
                @debug "Empty region -> Defaulting to leaf node"
                @debug "" min_val findall(new_nregs) findall(new_pregs) findall(reg_ids) branches
                jump_list[next] = 0 # pointing at root node -> leaf
                # Try to pick region that "disappeared", otherwise pick first region in parent
                reg_id  = isempty(fb_cands) ? findfirst(reg_ids) : findfirst(reg_ids .* .!new_regs_comp)
                hp_list[next] = fb_ids[reg_id]
            elseif length(fb_cands) > 1
                push!(U,(new_regs,branches ∪ [(hp_id,hp_sign)], next))
            else
                jump_list[next] = 0 # pointing at root node -> leaf
                hp_list[next] = first(fb_cands)
            end
        end
    end

    # Remove superfluous HPs
    hps,hp_list= remove_redundant_hps(jump_list,hp_list,hps)

    # Denormalize 
    hps = denormalize(hps,sol.scaling,sol.translation;hps=true)
    fbs = [denormalize(f,sol.scaling,sol.translation) for f in fbs]

    # Extract dual variables
    fbs_dual = dual ?  get_duals(CRs,sol) : Matrix{Float64}[] # XXX sparse representation makes C-code messy...

    # Cleanup
    DAQPBase.free_c_workspace(ws.p)

    zlims = clipping ? sol.problem.out_lims : zeros(0,2)
    return BinarySearchTree(hps,fbs,hp_list,jump_list,depth, fbs_dual, zlims)
end

function evaluate(bst::BinarySearchTree,θ)
    id =  1 # start in root note
    next_id = id+bst.jump_list[id]
    while next_id != id 
        hid = bst.hp_list[id]  
        if bst.halfplanes[1:end-1,hid]'*θ  ≤ bst.halfplanes[end,hid] 
            id = next_id+1 
        else
            id = next_id
        end
        next_id = id+bst.jump_list[id]
    end
    fid = bst.hp_list[id]  
    z = bst.feedbacks[fid]'*[θ;1]
    if(!isempty(bst.clipping))
        z = clamp.(z,bst.clipping[:,1],bst.clipping[:,2])
    end
    return z
end

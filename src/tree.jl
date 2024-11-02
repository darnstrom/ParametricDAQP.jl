struct BinarySearchTree
    halfplanes::Matrix{Float64}
    feedbacks::Vector{Matrix{Float64}}
    hp_list::Vector{Int}
    jump_list::Vector{Int}
end

function isnonempty(A,b)
    th,_,exitflag,_ = DAQP.quadprog(zeros(0,0),zeros(0),A,b,A_rowmaj = true)
    return exitflag == 1
end
function get_halfplanes(CRs)
    nreg = length(CRs)
    nreg == 0 && return nothing
    nth = size(CRs[1].Ath,1)
    hps = zeros(nth+1,0)
    reg2hp = [[] for _ in 1:nreg]

    for (reg_id,cr) in enumerate(CRs)
        for (i,a) = enumerate(eachcol(cr.Ath))
            # Disregard box bounds
            cr.bth[i] == 1 && sum(==(0),a) == 1 && continue
            asign = sign(a[findfirst(!=(0),a)])
            hcand = asign*[a;cr.bth[i]]
            # Check if hcand already exists in hps
            new_hp = true
            for (j,h) in  enumerate(eachcol(hps))
                if(all(hcand.≈h))
                    push!(reg2hp[reg_id],(j,asign))
                    new_hp = false
                    break;
                end
            end
            if new_hp
                hps = [hps hcand]
                push!(reg2hp[reg_id],(size(hps,2),asign))
            end
        end
    end
    return hps,reg2hp
end

function get_feedbacks(CRs; tol=1e-6)
    isempty(CRs) && return nothing
    nth = size(CRs[1].Ath,1)
    X,ids = [], []
    for cr in CRs 
        id = 0 
        for (i,x) in enumerate(X)
            if norm(cr.x-x) < tol
                id = i
                push!(ids,id)
                break
            end
        end
        if id == 0 # No duplicate
            push!(X,cr.x)
            push!(ids,length(X))
        end
    end
    return X,ids
end

# TODO: Can be cut in half by using points in CR 
function classify_regions(CRs,hps, reg2hp; reg_ids = nothing, hp_ids = nothing, branches = nothing)
    isnothing(reg_ids) && (reg_ids = 1:length(CRs))
    isnothing(hp_ids) && (hp_ids = 1:size(hps,2))
    nth, nh = size(hps,1)-1, length(hp_ids)
    
    nregs = [Set{Int}() for _ in 1:nh]
    pregs = [Set{Int}() for _ in 1:nh]


    Ath0,bth0 = zeros(nth,1),zeros(1)
    if !isnothing(branches)
        for (hid, hsign) in branches
            Ath0 = [Ath0 hsign*hps[1:nth,hid]]
            bth0 = [bth0; hsign*hps[end,hid]]
        end
    end

    for i in reg_ids 
        Ath_test = [Ath0 CRs[i].Ath]
        bth_test = [bth0;CRs[i].bth]
        for (j,hj) in enumerate(hp_ids)
            # First check if the hp is a facet of the region
            id = findfirst(x->first(x)==hj,reg2hp[i])
            if !isnothing(id) # hp is a facet of the region
                if last(reg2hp[i][id]) == 1
                    push!(pregs[j],i)
                else
                    push!(nregs[j],i)
                end
                continue
            end

            # Negative
            Ath_test[:,1] = -hps[1:nth,hj]
            bth_test[1] = -hps[end,hj]-(1e-6+1e-8)
            isnonempty(Ath_test,bth_test) && push!(nregs[j],i)
            # Positive
            Ath_test[:,1] = hps[1:nth,hj]
            bth_test[1] = hps[end,hj]-(1e-6+1e-8)
            isnonempty(Ath_test,bth_test) && push!(pregs[j],i)
        end
    end
    return nregs,pregs
end

function build_tree(CRs)
    hps,reg2hp = get_halfplanes(CRs)
    fbs, fb_ids = get_feedbacks(CRs)
    nh = size(hps,2)
    nregs,pregs = classify_regions(CRs,hps,reg2hp)
    hp_list, jump_list = Int[0],Int[0]

    N0 = (Set{Int}(1:length(CRs)),[],1)
    U = [N0]
    get_fbid = s->Set{Int}(fb_ids[collect(s)])
    split_objective = x-> max(length.(x)...) 
    while !isempty(U)
        reg_ids, branches, self_id = pop!(U) 
        # First use heuristic to find candidates 
        splits = [(reg_ids ∩ nregs[i], reg_ids ∩ pregs[i]) for i in 1:nh]
        vals = [split_objective(s) for s in splits]
        min_val = minimum(vals)
        hp_ids = findall(==(min_val),vals)
        if length(branches) > 0 # Not in root node -> have to compute the actual split
            splits = tuple.(classify_regions(CRs,hps,reg2hp;reg_ids = reg_ids, hp_ids = hp_ids, branches=branches)...)
            vals =[split_objective(s) for s in splits]
            min_val,min_id = findmin(vals)
            hp_id = hp_ids[min_id] # TODO add tie-breaker...
            new_nregs,new_pregs = splits[min_id]
        else
            hp_id = hp_ids[1] 
            new_nregs, new_pregs = splits[hp_id]
        end

        next_id = length(hp_list)+1
        hp_list[self_id] = hp_id
        jump_list[self_id] = next_id
        
        # Make room for the two new nodes that are spawned 
        push!(hp_list,0,0)
        push!(jump_list,0,0) 

        # Spawn new nodes
        fb_cands = get_fbid(new_nregs) 
        if length(fb_cands)!=1 
            push!(U,(new_nregs,branches ∪ [(hp_id,-1)], next_id))
        else
            jump_list[next_id] = 1 # pointing at root node -> leaf
            hp_list[next_id] = first(fb_cands) 
        end

        fb_cands = get_fbid(new_pregs) 
        if length(fb_cands)!=1 
            push!(U,(new_pregs,branches ∪ [(hp_id,1)], next_id+1))
        else
            jump_list[next_id+1] = 1 # pointing at root node -> leaf
            hp_list[next_id+1] = first(fb_cands) 
        end
    end
    return BinarySearchTree(hps,fbs,hp_list,jump_list)
end

function evaluate(bst::BinarySearchTree,θ)
    id =  1 # start in root note
    next_id = bst.jump_list[id]
    while next_id != 1
        hid = bst.hp_list[id]  
        if bst.halfplanes[1:end-1,hid]'*θ  ≤ bst.halfplanes[end,hid] 
            id = next_id+1 
        else
            id = next_id
        end
        next_id = bst.jump_list[id]
    end
    fid = bst.hp_list[id]  
    return bst.feedbacks[fid]'*[θ;1]
end

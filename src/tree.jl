struct BinarySearchTree
    halfplanes::Matrix{Float64}
    feedbacks::Vector{Matrix{Float64}}
    hp_list::Vector{Int}
    jump_list::Vector{Int}
    depth::Int
    duals::Vector{Matrix{Float64}}
end

function isnonempty(A,b;daqp_settings=nothing)
    d = DAQP.Model();
    DAQP.settings(d,Dict(:fval_bound=>size(A,1)-1,:sing_tol=>1e-9)) # Cannot be outside box
    !isnothing(daqp_settings) && DAQP.settings(d,daqp_settings)
    DAQP.setup(d,zeros(0,0),zeros(0),A,b,A_rowmaj=true);
    x,fval,exitflag,info = DAQP.solve(d);
    # TODO add guard to do this check
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
            cr.bth[i] == 1 && sum(!=(0),a) == 1 && continue
            nz_id = findfirst(!=(0),a)
            isnothing(nz_id) && continue
            asign = sign(a[nz_id])
            hcand = asign*[a;cr.bth[i]]
            # Check if hcand already exists in hps
            new_hp = true
            for (j,h) in  enumerate(eachcol(hps))
                if(all(isapprox.(hcand,h,atol=1e-5,rtol=1e-5)))
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

function get_feedbacks(CRs; tol=1e-5)
    isempty(CRs) && return nothing
    nth = size(CRs[1].Ath,1)
    Z,ids = [], []
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
function classify_regions(CRs,hps, reg2hp; reg_ids = nothing, hp_ids = nothing, branches = nothing, daqp_settings=nothing)
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
        bth_test = [bth0;CRs[i].bth].-1e-6
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
            bth_test[1] = -hps[end,hj]-1e-6
            isnonempty(Ath_test,bth_test;daqp_settings) && push!(nregs[j],i)
            # Positive
            Ath_test[:,1] = hps[1:nth,hj]
            bth_test[1] = hps[end,hj]-1e-6
            isnonempty(Ath_test,bth_test;daqp_settings) && push!(pregs[j],i)
        end
    end
    return nregs,pregs
end

function build_tree(sol::Solution; daqp_settings = nothing, verbose=1, max_reals=1e12, dual=false)
    if sol.status != :Solved
        verbose > 0 && @warn "Cannot build binary search tree. Solution status: $(sol.status)"
        return nothing
    end
    sol.status != :Solved && return nothing
    verbose > 0 && @info "Building binary search tree" 
    hps,reg2hp = get_halfplanes(sol.CRs)
    fbs, fb_ids = get_feedbacks(sol.CRs)

    # Approximate number of real numbers
    n_reals = length(hps) + length(fbs)*length(fbs[1])
    dual && (n_reals += length(sol.CRs)*(sol.problem.n_theta+1)*length(sol.problem.norm_factors))
    if n_reals > max_reals
        verbose > 0 && @warn "Memory limit for real numbers (n_reals = $n_reals, max_reals = $max_reals) reached"
        return nothing
    end

    nh = size(hps,2)
    nregs,pregs = classify_regions(sol.CRs,hps,reg2hp;daqp_settings)
    hp_list, jump_list = Int[0],Int[0]

    N0 = (Set{Int}(1:length(sol.CRs)),[],1)
    U = [N0]
    get_fbid = s->Set{Int}(fb_ids[collect(s)])
    split_objective = dual ? x->max(length.(x)...) : x-> max(length.(get_fbid.(x))...)

    depth = 0
    while !isempty(U)
        reg_ids, branches, self_id = pop!(U)
        depth = max(depth,length(branches))

        hp_ids = reduce(∪,Set(first.(reg2hp[i])) for i in reg_ids);
        hp_ids = collect(setdiff!(hp_ids,first(b) for b in branches))

        # First use heuristic to find candidates 
        splits = [(reg_ids ∩ nregs[i], reg_ids ∩ pregs[i]) for i in hp_ids]
        vals = [split_objective(s) for s in splits]
        min_val = minimum(vals)
        min_ids = findall(==(min_val),vals)
        hp_ids = hp_ids[min_ids]
        if length(branches) > 0 && min_val > 1# Compute the actual split
            splits = tuple.(classify_regions(sol.CRs,hps,reg2hp;reg_ids,hp_ids,branches,daqp_settings)...)
            vals =[split_objective(s) for s in splits]
            min_val = minimum(vals)
            min_ids = findall(==(min_val),vals)
            # Among the minimum values, do new split to maximum regions that split
            if(!dual)
                vals =[max(length.(s)...) for s in splits[min_ids]]
                min_val_second,min_id = findmin(vals)

                min_ids = min_ids[findall(==(min_val_second),vals)]
                vals =[min(length.(get_fbid.(s))...) for s in splits[min_ids]]
                min_val_third,min_id = findmin(vals)
            else
                min_id = 1; # Pick first
            end

            hp_id = hp_ids[min_ids[min_id]]
            new_nregs,new_pregs = splits[min_ids[min_id]]
        else
            hp_id = hp_ids[1] 
            new_nregs, new_pregs = splits[min_ids[1]]
        end

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
                @warn "Empty region -> Defaulting to leaf node" min_val new_nregs new_pregs reg_ids branches
                jump_list[next] = 0 # pointing at root node -> leaf
                # Try to pick region that "disappeared", otherwise pick first region in parent
                reg_id  = isempty(fb_cands) ? first(reg_ids) : first(setdiff(reg_ids,new_regs_comp))
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
    hp_rows = findall(jump_list.!== 0)
    new_hp_ids = sort(unique(hp_list[hp_rows]))
    hp_old2new = Dict(id => k for (k,id) in enumerate(new_hp_ids))
    new_hps = hps[:,new_hp_ids]

    new_hp_list = copy(hp_list)
    for i in hp_rows
        new_hp_list[i] = hp_old2new[new_hp_list[i]]
    end

    # Denormalize 
    hps = denormalize(new_hps,sol.scaling,sol.translation;hps=true)
    fbs = [denormalize(f,sol.scaling,sol.translation) for f in fbs]

    fbs_dual = Matrix{Float64}[]
    if dual # XXX would prefer sparse representation, but makes C-code a lot messier...
        for cr in sol.CRs
            lam = zeros(size(cr.lam,1),length(sol.problem.norm_factors))
            lam[:,cr.AS] = cr.lam
            push!(fbs_dual,lam)
        end
        fbs_dual = [denormalize(f,sol.scaling,sol.translation) for f in fbs_dual]
    end

    return BinarySearchTree(hps,fbs,new_hp_list,jump_list,depth, fbs_dual)
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
    return bst.feedbacks[fid]'*[θ;1]
end

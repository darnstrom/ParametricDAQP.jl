const _loaded_worker_ids = Set{Int}()
const _explicit_context = Ref{Any}(nothing)
const _explicit_workspace = Ref{Any}(nothing)
const _classify_context = Ref{Any}(nothing)
const _classify_workspace = Ref{Any}(nothing)

const _MIN_EXPLICIT_BATCH_PER_WORKER = 4
const _MIN_CLASSIFY_WORK_PER_WORKER = 64

function _distributed_workers()
    Distributed.myid() != 1 && return Int[]
    pids = filter(!=(Distributed.myid()), Distributed.workers())
    isempty(pids) && return pids

    stale_ids = setdiff(_loaded_worker_ids, Set(pids))
    foreach(pid -> delete!(_loaded_worker_ids, pid), stale_ids)
    for pid in pids
        pid in _loaded_worker_ids && continue
        fetch(Distributed.remotecall_eval(Main, pid, :(using ParametricDAQP)))
        push!(_loaded_worker_ids, pid)
    end
    return pids
end

function _split_batches(items, nbatches)
    isempty(items) && return Vector{typeof(items)}()
    nbatches = max(1, min(nbatches, length(items)))
    batch_size = cld(length(items), nbatches)
    return [items[i:min(i + batch_size - 1, end)] for i in 1:batch_size:length(items)]
end

function _should_parallelize_explicit(npending, nworkers)
    return nworkers > 0 && npending >= nworkers * _MIN_EXPLICIT_BATCH_PER_WORKER
end

function _should_parallelize_classification(nreg, nhp, nworkers)
    return nworkers > 0 && nreg > 1 && nhp > 0 && nreg * nhp >= nworkers * _MIN_CLASSIFY_WORK_PER_WORKER
end

function _clear_explicit_context!()
    if !isnothing(_explicit_workspace[])
        DAQPBase.free_c_workspace(_explicit_workspace[].DAQP_workspace)
        _explicit_workspace[] = nothing
    end
    _explicit_context[] = nothing
    return nothing
end

function _set_explicit_context!(prob, Θ, opts)
    _clear_explicit_context!()
    _explicit_context[] = (prob=prob, Θ=Θ, opts=opts)
    return nothing
end

function _get_explicit_workspace()
    ws = _explicit_workspace[]
    if isnothing(ws)
        ctx = _explicit_context[]
        ws = setup_workspace(ctx.Θ, length(ctx.prob.norm_factors); opts=ctx.opts)
        ws.is_equality[ctx.prob.eq_ids] .= true
        _explicit_workspace[] = ws
    end
    return ws
end

function _process_explicit_batch(batch)
    isempty(batch) && return (results=NamedTuple[], nLPs=0)
    ctx = _explicit_context[]
    ws = _get_explicit_workspace()
    nLPs0 = ws.nLPs
    T = eltype(batch)
    results = Vector{NamedTuple{(:as, :region, :up, :down), Tuple{T, Union{Nothing, CriticalRegion}, Bool, Bool}}}(undef, length(batch))

    for (k, as) in enumerate(batch)
        region, up, down = isoptimal(as, ws, ctx.prob, ctx.opts)
        results[k] = (as=as, region=region, up=up, down=down)
    end
    return (results=results, nLPs=ws.nLPs - nLPs0)
end

function _initialize_explicit_workers!(prob, Θ, opts)
    pids = _distributed_workers()
    for pid in pids
        Distributed.remotecall_fetch(ParametricDAQP._set_explicit_context!, pid, prob, Θ, opts)
    end
    return pids
end

function _clear_explicit_workers!(pids)
    for pid in pids
        Distributed.remotecall_fetch(ParametricDAQP._clear_explicit_context!, pid)
    end
    return nothing
end

function _parallel_process_explicit_batches(pending, pids)
    batches = _split_batches(pending, length(pids))
    pool = Distributed.CachingPool(pids)
    return Distributed.pmap(ParametricDAQP._process_explicit_batch, pool, batches)
end

function _clear_classify_context!()
    if !isnothing(_classify_workspace[])
        DAQPBase.free_c_workspace(_classify_workspace[].p)
        _classify_workspace[] = nothing
    end
    _classify_context[] = nothing
    return nothing
end

function _set_classify_context!(CRs, hps, reg2hp, n_theta, max_region_constraints)
    _clear_classify_context!()
    _classify_context[] = (
        CRs=CRs,
        hps=hps,
        reg2hp=reg2hp,
        n_theta=n_theta,
        max_region_constraints=max_region_constraints
    )
    return nothing
end

function _get_classify_workspace(n_branches)
    ctx = _classify_context[]
    ws = _classify_workspace[]
    nc_max = max(250, 1 + n_branches + ctx.max_region_constraints)
    if isnothing(ws) || length(ws.b) < nc_max
        !isnothing(ws) && DAQPBase.free_c_workspace(ws.p)
        ws = setup_daqp_workspace(ctx.n_theta, nc_max)
        _classify_workspace[] = ws
    end
    return ws
end

function _process_classify_chunk(args)
    reg_ids, hp_ids, branches = args
    ctx = _classify_context[]
    ws = _get_classify_workspace(isnothing(branches) ? 0 : length(branches))
    nh = length(hp_ids)
    nregs = falses(nh, length(reg_ids))
    pregs = falses(nh, length(reg_ids))
    _classify_regions_chunk!(nregs, pregs, ctx.CRs, ctx.hps, ctx.reg2hp, ws, reg_ids, hp_ids, branches; local_columns=true)
    return (reg_ids=reg_ids, nregs=nregs, pregs=pregs)
end

function _initialize_classify_workers!(CRs, hps, reg2hp, n_theta, max_region_constraints)
    pids = _distributed_workers()
    for pid in pids
        Distributed.remotecall_fetch(
            ParametricDAQP._set_classify_context!,
            pid,
            CRs,
            hps,
            reg2hp,
            n_theta,
            max_region_constraints
        )
    end
    return pids
end

function _clear_classify_workers!(pids)
    for pid in pids
        Distributed.remotecall_fetch(ParametricDAQP._clear_classify_context!, pid)
    end
    return nothing
end

function _parallel_classify_regions(CRs, hps, reg2hp, reg_ids, hp_ids, branches, n_theta, nR, ws)
    pids = _distributed_workers()
    batches = _split_batches(reg_ids, length(pids) + 1)
    isempty(batches) && return [falses(nR) for _ in hp_ids], [falses(nR) for _ in hp_ids]

    payload = [(batch, hp_ids, branches) for batch in batches]
    chunk_results = Vector{Any}(undef, length(payload))
    remote_jobs = min(length(pids), max(0, length(payload) - 1))
    futures = Vector{Any}(undef, remote_jobs)
    for i in 1:remote_jobs
        futures[i] = Distributed.remotecall(ParametricDAQP._process_classify_chunk, pids[i], payload[i])
    end

    local_payload = payload[end]
    local_nregs = falses(length(hp_ids), length(local_payload[1]))
    local_pregs = falses(length(hp_ids), length(local_payload[1]))
    _classify_regions_chunk!(local_nregs, local_pregs, CRs, hps, reg2hp, ws, local_payload[1], hp_ids, branches; local_columns=true)
    chunk_results[end] = (reg_ids=local_payload[1], nregs=local_nregs, pregs=local_pregs)
    for i in 1:remote_jobs
        chunk_results[i] = fetch(futures[i])
    end

    nregs_mat = falses(length(hp_ids), nR)
    pregs_mat = falses(length(hp_ids), nR)
    for chunk in chunk_results
        for (local_col, reg_id) in enumerate(chunk.reg_ids)
            @views nregs_mat[:, reg_id] .= chunk.nregs[:, local_col]
            @views pregs_mat[:, reg_id] .= chunk.pregs[:, local_col]
        end
    end
    return _rows_to_bitvectors(nregs_mat), _rows_to_bitvectors(pregs_mat)
end

## Plots.jl
using RecipesBase
@recipe function f(r::CriticalRegion)
    st --> :mesh3d 
    nth = size(r.Ath,1)
    legend --> false
    vs=PolyDAQP.vrep_2d(minrep(slice(r.Ath,r.bth,collect(3:nth))...)...)
    nv = length(vs)
    x,y = first.(vs),last.(vs)
    z = [r.xTH[:,1]'*[v;zeros(nth-2)] + r.xC[1] for v in vs]
    connections--> (zeros(Int, nv-2),collect(1:nv-2),collect(2:nv-1))
    return x,y,z 
end

@recipe function f(rs::Vector{CriticalRegion};uid=0, 
        fix_ids = nothing, fix_vals=nothing)
    isempty(rs) && error("Cannot plot empty collection")
    nth = size(rs[1].Ath,1)
    ids = isnothing(fix_ids) ? collect(3:nth) : fix_ids
    values = isnothing(fix_vals) ? zeros(nth-2) : fix_vals 
    free_ids = setdiff(1:nth,ids) 
    if(uid ==0) # Don't plot feedback law, just the poly collection
        @series [Polyhedron(slice(r.Ath,r.bth,ids;values)...) for r in rs]
    else
        for r in rs
            vs=PolyDAQP.vrep_2d(minrep(slice(r.Ath,r.bth,ids;values)...)...)
            nv = length(vs)
            nv < 2 && continue
            x,y = first.(vs), last.(vs)
            c = r.xTH[ids,uid]'*values + r.xC[uid] 
            z = [r.xTH[free_ids,uid]'*v + c for v in vs]
            @series begin
                st --> :mesh3d 
                legend --> false
                connections--> (zeros(Int, nv-2),collect(1:nv-2),collect(2:nv-1))
                extra_kwargs --> Dict(:subplot=>Dict("faceted color" => "none"))
                (x,y,z)
            end
        end
    end
end

## PGFPlotsX
function pplot(rs::Vector{CriticalRegion};uid=0, fix_ids = nothing, fix_vals=nothing,opts=Dict{Symbol,Any}())
    isempty(rs) && error("Cannot plot empty collection")
    nth = size(rs[1].Ath,1)
    ids = isnothing(fix_ids) ? collect(3:nth) : fix_ids
    values = isnothing(fix_vals) ? zeros(nth-2) : fix_vals
    free_ids = setdiff(1:nth,ids)

    ps = PolyDAQP.Polyhedron[]
    fs = Function[]
    for r in rs
        (!isempty)
        p = Polyhedron(slice(r.Ath,r.bth,ids;values)...)
        isempty(p) && continue
        push!(ps,p)
        uid == 0 && continue
        c = r.xTH[ids,uid]'*values + r.xC[uid]
        push!(fs,v->c+r.xTH[free_ids,uid]'*v)
    end
    PolyDAQP.pplot(ps;fs,opts)
end

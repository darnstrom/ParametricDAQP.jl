## Plots.jl
using RecipesBase
@recipe function f(r::CriticalRegion)
    st --> :mesh3d 
    nth = size(r.Ath,1)
    legend --> false
    vs=PolyDAQP.vrep_2d(minrep(slice(r.Ath,r.bth,collect(3:nth))...)...)
    nv = length(vs)
    x,y = first.(vs),last.(vs)
    z = [r.x[:,1]'*[v;zeros(nth-2);1] for v in vs]
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
            c = r.x[ids,uid]'*values + r.x[end,uid]
            z = [r.x[free_ids,uid]'*v + c for v in vs]
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
function plot_regions(sol::Solution;fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    pplot(get_critical_regions(sol);out_id=0,fix_ids,fix_vals,opts)
end
function plot_solution(sol::Solution;out_id=1,fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    pplot(get_critical_regions(sol);out_id,fix_ids,fix_vals,opts)
end
function pplot(rs::Vector{CriticalRegion};out_id=0, fix_ids = nothing, fix_vals=nothing,opts=Dict{Symbol,Any}())
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
        out_id == 0 && continue
        c = r.x[ids,out_id]'*values + r.x[end,out_id]
        push!(fs,v->c+r.x[free_ids,out_id]'*v[1:2])
    end
    # Some plotting
    lopts = Dict(
                 :xlabel=>"\\large\$\\theta_{"*string(free_ids[1])*"}\$",
                 :ylabel=>"\\large\$\\theta_{"*string(free_ids[2])*"}\$",
                )
    if out_id != 0
        push!(lopts,:zlabel=>"\\large\$x^*_{"*string(out_id[1])*"}\$")

        lopts = merge(Dict(:view=>(45,45),),lopts)
    else
        push!(lopts,:ylabel_style=> "{yshift={-5pt}}")
    end

    opts = merge(lopts,opts)
    PolyDAQP.pplot(ps;fs,opts)
end

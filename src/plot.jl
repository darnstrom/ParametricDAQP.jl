## Plots.jl
using RecipesBase
@recipe function f(r::CriticalRegion; z_id =0, free_ids = zeros(0), fix_ids = zeros(0), fix_vals=zeros(0))
    plotattributes[:CR_attr] = (z_id ,free_ids,fix_ids,fix_vals)
    return [r]
end

@recipe function f(rs::Vector{CriticalRegion};z_id =0, free_ids=zeros(0), fix_ids = zeros(0), fix_vals=zeros(0))
    isempty(rs) && error("Cannot plot empty collection")
    nth = size(rs[1].Ath,1)

    if haskey(plotattributes, :CR_attr)
        z_id ,free_ids,fix_ids,fix_vals = pop!(plotattributes,:CR_attr)
    end

    if isempty(free_ids)
        ids = isempty(fix_ids) ? collect(3:nth) : fix_ids
        values = isempty(fix_vals) ? zeros(nth-2) : fix_vals
        free_ids = setdiff(1:nth,ids)
    elseif length(free_ids) != 2
        error("The number of parameters to plot needs to be 2, not $(length(free_ids))")
    else
        ids = setdiff(1:nth,free_ids)
        values = isempty(fix_vals) ? zeros(nth-2) : fix_vals
    end
    free_ids = sort(free_ids)
    xlabel --> "\\theta [$(free_ids[1])]"
    ylabel --> "\\theta [$(free_ids[2])]"
    ps = [Polyhedron(slice(r.Ath,r.bth,ids;values)...) for r in rs]
    if(z_id  ==0) # Don't plot feedback law, just the poly collection
        title --> "Critical regions"
        return ps
    else
        title --> "Explicit solution"
        zlabel --> "z [$z_id]"
        fs = [(v->r.z[free_ids,z_id ]'*v+r.z[ids,z_id ]'*values + r.z[end,z_id ]) for r in rs]
        return collect(zip(ps,fs))
    end
end

@recipe function f(sol::Solution; z_id  = 0, free_ids=zeros(0), fix_ids = zeros(0), fix_vals=zeros(0))
    plotattributes[:CR_attr] = (z_id ,free_ids,fix_ids,fix_vals)
    return get_critical_regions(sol)
end

## PGFPlotsX
function plot_regions(sol::Solution;fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    pplot(get_critical_regions(sol);z_id=0,fix_ids,fix_vals,opts)
end
function plot_solution(sol::Solution;z_id=1,fix_ids=nothing,fix_vals=nothing,opts=Dict{Symbol,Any}())
    pplot(get_critical_regions(sol);z_id,fix_ids,fix_vals,opts)
end
function pplot(rs::Vector{CriticalRegion};z_id=0, fix_ids = nothing, fix_vals=nothing,opts=Dict{Symbol,Any}())
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
        z_id == 0 && continue
        c = r.z[ids,z_id]'*values + r.z[end,z_id]
        push!(fs,v->c+r.z[free_ids,z_id]'*v[1:2])
    end
    # Some plotting
    lopts = Dict(
                 :xlabel=>"\\large\$\\theta_{"*string(free_ids[1])*"}\$",
                 :ylabel=>"\\large\$\\theta_{"*string(free_ids[2])*"}\$",
                )
    if z_id != 0
        push!(lopts,:zlabel=>"\\large\$z^*_{"*string(z_id[1])*"}\$")

        lopts = merge(Dict(:view=>(45,45),),lopts)
    else
        push!(lopts,:ylabel_style=> "{yshift={-5pt}}")
    end

    opts = merge(lopts,opts)
    PolyDAQP.pplot(ps;fs,opts)
end

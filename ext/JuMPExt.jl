module JuMPExt

using JuMP 
using ParametricDAQP

function ParametricDAQP.mpsolve(model::JuMP.Model,vars;
        out_vars=nothing,opts=nothing, eliminate_equalities=true)
    ParametricDAQP.mpsolve(model.moi_backend,vars;out_vars,opts,eliminate_equalities)
end
end

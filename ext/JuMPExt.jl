module JuMPExt

using JuMP 
using ParametricDAQP

function ParametricDAQP.mpsolve(model::JuMP.Model,vars;
        out_vars=nothing,opts=nothing, eliminate_equalities=true)
    ParametricDAQP.mpsolve(model.moi_backend,vars;out_vars,opts,eliminate_equalities)
end

function ParametricDAQP.get_mpp(model::JuMP.Model,vars;out_vars=nothing, eliminate_equalities=true)
    ParametricDAQP.get_mpp(model.moi_backend,vars;out_vars,eliminate_equalities)
end

function ParametricDAQP.codegen_implicit(model::JuMP.Model, vars; out_vars=nothing, eliminate_equalities=false, fname="pdaqp_workspace", dir="codegen", opt_settings=nothing, src=true, float_type="double",warm_start=false)
    ParametricDAQP.codegen_implicit(model.moi_backend, vars; out_vars, eliminate_equalities, fname, dir, opt_settings, src, float_type,warm_start)
end

end

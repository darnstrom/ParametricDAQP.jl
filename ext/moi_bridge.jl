using DAQPBase 
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

const Interval = MOI.Interval{Cdouble}

const Affine = MOI.ScalarAffineFunction{Cdouble}
const Quadratic = MOI.ScalarQuadraticFunction{Cdouble}
const VectorAffine = MOI.VectorAffineFunction{Cdouble}

const Interval = MOI.Interval{Cdouble}
const LessThan = MOI.LessThan{Cdouble}
const GreaterThan = MOI.GreaterThan{Cdouble}
const EqualTo = MOI.EqualTo{Cdouble}
const SupportedSets= Union{GreaterThan, LessThan, EqualTo, Interval}


mutable struct ParametricOptimizer <: MOI.AbstractOptimizer
    sense::MOI.OptimizationSense
    objconstant::Cdouble
    rows::Dict{Int, Int}
    is_empty::Bool

    function ParametricOptimizer()
        sense = MOI.MIN_SENSE
        objconstant = 0.0 
        rows = Dict{Int, Int}()
        return new(sense,objconstant,rows,true)
    end
end

MOI.is_empty(optimizer::ParametricOptimizer) = optimizer.is_empty 
function MOI.empty!(optimizer::ParametricOptimizer)
    optimizer.is_empty = true
    optimizer.sense = MOI.MIN_SENSE
    optimizer.objconstant = 0.0 
    optimizer.rows = Dict{Int, Int}()
end

## Supported constraint types

MOI.supports_constraint(
    ::ParametricOptimizer,
    ::Type{<:MOI.ScalarAffineFunction},
    ::Type{<:SupportedSets}
) = true

MOI.supports_constraint(
    ::ParametricOptimizer,
    ::Type{<:MOI.VariableIndex},
    ::Type{<:SupportedSets}
) = true

## Supported objective functions
MOI.supports(
    ::ParametricOptimizer,
    ::MOI.ObjectiveFunction{<:Union{
       MOI.ScalarAffineFunction,
       MOI.ScalarQuadraticFunction,
    }}
) = true

## Helper functions 
function check_attributes(dest::ParametricOptimizer, src)
    #allowable model attributes
    for attr in MOI.get(src, MOI.ListOfModelAttributesSet())
        if attr == MOI.Name()           ||
            attr == MOI.ObjectiveSense() ||
            attr isa MOI.ObjectiveFunction
            continue
        end
        throw(MOI.UnsupportedAttribute(attr))
    end

    #allowable variable attributes
    for attr in MOI.get(src, MOI.ListOfVariableAttributesSet())
        if attr == MOI.VariableName() || attr == MOI.VariablePrimalStart()
            continue
        end
        throw(MOI.UnsupportedAttribute(attr))
    end

    #allowable constraint types and attributes
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        if !MOI.supports_constraint(dest, F, S)
            throw(MOI.UnsupportedConstraint{F, S}())
        end
        for attr in MOI.get(src, MOI.ListOfConstraintAttributesSet{F, S}())
            if attr == MOI.ConstraintName()
                continue
            end
            throw(MOI.UnsupportedAttribute(attr))
        end
    end

    return nothing
end

# Set up index map from `src` variables and constraints to `dest` variables and constraints.
function MOIU.IndexMap(dest::ParametricOptimizer, src::MOI.ModelLike)

    idxmap = MOIU.IndexMap()

    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = MOI.VariableIndex(i)
    end
    i = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        MOI.supports_constraint(dest, F, S) || throw(MOI.UnsupportedConstraint{F, S}())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci in cis_src
            i += 1
            idxmap[ci] = MOI.ConstraintIndex{F, S}(i)
        end
    end

    return idxmap
end


function assign_constraint_rows!(
    dest::ParametricOptimizer,
    idxmap::MOIU.IndexMap,
    src::MOI.ModelLike
)

    startrow = 1+MOI.get(src, MOI.NumberOfVariables())
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci_src in cis_src
            ci_dest = idxmap[ci_src]
            if(F!=MOI.VariableIndex)
                dest.rows[ci_dest.value] = startrow 
                startrow += 1
            else
                var_id = idxmap[MOI.get(src, MOI.ConstraintFunction(), ci_src)].value
                dest.rows[ci_dest.value] = var_id
            end
        end
    end

    return nothing
end

## Extract Constraints
# blower[1:n]     ≤  x ≤  bupper[1:n] 
# blower[n+1:end] ≤ Ax ≤  bupper[n+1:end]
function process_constraints(
    dest::ParametricOptimizer,
    src::MOI.ModelLike,
    idxmap
)

    rows = dest.rows
    n = MOI.get(src, MOI.NumberOfVariables())
    m = (isempty(rows)) ? n : maximum(values(rows))

    A = zeros(Cdouble,n,m-n)
    bupper = Vector{Cdouble}(undef, m)
    blower = Vector{Cdouble}(undef, m)
    offset = Vector{Cdouble}(undef, m)
    sense  = Vector{Cint}(undef,m)

    bupper[1:n].=1e30;
    blower[1:n].=-1e30;
    offset[1:n].=0;
    sense[1:n].= DAQPBase.IMMUTABLE


    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        process_constraints!(A, bupper, blower, sense, offset,
                            src, idxmap, rows, 
                            F, S)
    end
    bupper .-= offset 
    blower .-= offset 

    return (A, bupper, blower, sense)

end

function process_constraints!(
    A::Matrix{Cdouble},
    bupper::Vector{Cdouble},
    blower::Vector{Cdouble},
    sense::Vector{Cint},
    offset::Vector{Cdouble},
    src::MOI.ModelLike,
    idxmap,
    rows::Dict{Int,Int},
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
)
    n = MOI.get(src, MOI.NumberOfVariables())
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        row = rows[idxmap[ci].value]
        extract_offset(offset, row, f)
        extract_A(A, f, row.-n, idxmap)
        extract_b(bupper, blower, sense, row, f, s)
    end
    return
end

extract_offset(::Vector{Cdouble}, ::Int, ::MOI.VariableIndex) = nothing
extract_A(::Matrix{Cdouble},::MOI.VariableIndex,::Int, Any) = nothing

function extract_offset(offset::Vector{Cdouble}, row::Int, f::Affine)
    offset[row] = MOI.constant(f, Cdouble)
    return
end

function extract_A(A::Matrix{Cdouble}, f::Affine, row::Int, idxmap)
    for term in f.terms
        var = term.variable
        col = idxmap[var].value
        A[col,row] = term.coefficient # colmaj -> rowmaj
    end
end

function extract_b(
    bupper::Vector{Cdouble},
    blower::Vector{Cdouble},
    sense::Vector{Cint},
    row::Int,
    f::Affine,
    s::SupportedSets,
)
    extract_b(bupper,blower,sense, row, MOI.Interval(s))
    return
end

function extract_b(
    bu::Vector{Cdouble},
    bl::Vector{Cdouble},
    sense::Vector{Cint},
    row::Int,
    f::MOI.VariableIndex,
    s::SupportedSets,
)

    i = (s==MOI.ZeroOne()) ? MOI.Interval(0,1) : MOI.Interval(s)
    extract_b(bu,bl,sense, row, MOI.Interval(max(i.lower,bl[row]),min(i.upper,bu[row])))
    if(s==MOI.ZeroOne()) sense[row] = BINARY end # Mark binary constraints
end

function extract_b(
    bupper::Vector{Cdouble},
    blower::Vector{Cdouble},
    sense::Vector{Cint},
    row::Int,
    interval::Interval,
)
    bupper[row] = interval.upper
    blower[row] = interval.lower
    sense[row] = (interval.lower == interval.upper) ?  DAQPBase.EQUALITY : 0
    return
end

## Extract Objective
# Construct the objective minimize `0.5 x' H x + f' x + c
function process_objective(dest::ParametricOptimizer, src::MOI.ModelLike, idxmap)
    sense = dest.sense
    n = MOI.get(src, MOI.NumberOfVariables())

    if sense == MOI.FEASIBILITY_SENSE
        H,f,c = zeros(n,n),zeros(n),0.0
    else
        function_type = MOI.get(src, MOI.ObjectiveFunctionType())
        f = zeros(Cdouble,n)

        if function_type == Affine
            faffine = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction}())
            H = zeros(n,n)
            process_objective_linearterm!(f, faffine.terms, idxmap)
            c = faffine.constant

        elseif function_type == Quadratic 
            H = zeros(Cdouble,n,n)
            fquadratic = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction}())
            for term in fquadratic.quadratic_terms
                i,j = Int(idxmap[term.variable_1].value),Int(idxmap[term.variable_2].value)
                H[i,j] += term.coefficient
                if(i!=j)
                    H[j,i] += term.coefficient
                end
            end
            process_objective_linearterm!(f, fquadratic.affine_terms, idxmap)
            c = fquadratic.constant

        else
            throw(MOI.UnsupportedAttribute(MOI.ObjectiveFunction{function_type}()))
        end

        if sense == MOI.MAX_SENSE
            H,f,c = -H,-f,-c
        end

    end
    return (H, f, c)
end

function process_objective_linearterm!(
    f::Vector{Cdouble},
    terms::Vector{<:MOI.ScalarAffineTerm},
    idxmapfun::Function = identity
) 
    f .= 0
    for term in terms
        var = term.variable
        coeff = term.coefficient
        f[idxmapfun(var).value] += coeff
    end
    return nothing
end

function process_objective_linearterm!(
    f::Vector{Cdouble},
    terms::Vector{<:MOI.ScalarAffineTerm},
    idxmap::MOIU.IndexMap
)
    process_objective_linearterm!(f, terms, var -> idxmap[var])
end

function generate_mpQP(n,m,nth)
    M = randn(n,n)
    H = M*M'
    f = randn(n,1)
    F = randn(n,nth)
    A = randn(m,n)
    b = [rand(m,1);rand(m,1)]
    F0 = randn(n,nth); # The point F0*th will be primal feasible
    B =[A;-A]*(-F0);
    #bounds_table = [collect(m+1:2m);collect(1:m)]
    bounds_table = collect(1:2m)
    senses = zeros(Cint,2m)
    mpQP = (H=H,f=f,F=F,H_theta=zeros(0,0),
                A=[A;-A],b=b,B=B,bounds_table=bounds_table,senses=senses)

    P_theta = (A = zeros(nth,0), b=zeros(0), ub=ones(nth),lb=-ones(nth),F0=F0) 

    return mpQP,P_theta
end
function generate_mpVI(n, m, nth)
    H = randn(n, n)
    eig_min = minimum(eigvals(H + H'))
    if eig_min < 0
        H = H - 2 * eig_min * I(n)
    end
    f = randn(n)
    F = randn(n, nth)
    A = randn(m, n)
    b = rand(2 * m)
    F0 = randn(n, nth) # The point F0*th will be primal feasible
    B = [A; -A] * (-F0)

    mpVI = ParametricDAQP.MPVI(H, F, f, [A; -A], B, b)
    P_theta = (A=zeros(nth, 0), b=zeros(0), ub=ones(nth), lb=-ones(nth), F0=F0)

    return mpVI, P_theta
end

function pointlocation(th::Vector{Float64}, partition ;eps_gap=1e-6)
    contained_in= Int64[]
    for (ind,region) in enumerate(partition)
        violation = minimum(region.bth-region.Ath'*th)
        if(violation>=-eps_gap)
            push!(contained_in,ind)
        end
    end
    return contained_in
end
function evaluate_solution(sol::ParametricDAQP.Solution, θ::Vector{Float64}; eps_gap=1e-6)
    θ_normalized = (θ - sol.translation) .* sol.scaling
    all_CRs_indexes = pointlocation(θ_normalized, sol.CRs)
    if !isempty(all_CRs_indexes)
        CR = sol.CRs[all_CRs_indexes[1]]
        return CR.z' * [θ_normalized; 1]
    else
        return nothing
    end
end

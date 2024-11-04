using ParametricDAQP
using PolyDAQP
using Test
using LinearAlgebra
using DAQP

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

@testset "ParametricDAQP.jl" begin
    n,m,nth = 5,20,6
    tol = 1e-5

    opts = ParametricDAQP.Settings()
    opts.store_points=true
    opts.store_dual=true


    mpQP,Θ = generate_mpQP(n,m,nth)
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);

    # Test for 10000 random points 
    N = 10000
    ths = 2*rand(nth,N).-1;
    containment_inds = zeros(Int,N) 
    errs_primal,errs_dual = zeros(N),zeros(N)
    for n = 1:N
        θ = ths[:,n]
        inds = pointlocation(θ,sol.CRs);
        containment_inds[n]=length(inds)
        AS = sol.CRs[inds[1]].AS
        xsol = sol.CRs[inds[1]].x'*[θ;1]
        λsol = sol.CRs[inds[1]].lam'*[θ;1]
        f = mpQP.f[:,1]+mpQP.F*θ
        b = mpQP.b[:,1]+mpQP.B*θ
        xref,~,~,info= DAQP.quadprog(mpQP.H,f,mpQP.A,b,-1e30*ones(2m),mpQP.senses);
        errs_primal[n] = norm(xsol-xref) 
        errs_dual[n] = norm(λsol-info.λ[sol.CRs[inds[1]].AS]) 
    end
    @test all(errs_primal.< tol)
    @test all(errs_dual.< tol)
    @test ~any(containment_inds.==0) # No holes
    @test sum(containment_inds.>1) < 0.01*N # Not more than 1% overlap
end

# Example 2 in Ahmadi-Moshkenani 2018
@testset "Full-dim LICQ" begin 
    n, nth = 4,2
    H = diagm(ones(n)) 
    f = zeros(n,1)
    F = zeros(n,nth)
    A = [[1 0;
         -1 0;
         0 1;
         0 -1;] -ones(4,1) 0.5*ones(4,1)]
    b = -ones(4,1)
    B = [1.0 0; -1 0; 0 -1; 1 1];
    mpQP = (H=H,f=f,F=F,
            A=A,b=b,B=B)
    Θ = (ub=2*ones(nth),lb=-2*ones(nth)) 

    opts = ParametricDAQP.Settings()
    opts.verbose=1;
    opts.postcheck_rank=true
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
end

# Non facet to facet
@testset "Violated Facet-to-facet" begin 
    n, nth = 3,2
    H = diagm(ones(n)) 
    f = zeros(n,1)
    F = zeros(n,nth)
    A = [1 0 -1;
         -1 0 -1;
         0 1 -1;
         0 -1 -1;
         3/4 16/25 -1;
         -3/4 -16/25 -1]
    b = -ones(6);
    B = [1.0 0; -1 0; 0 -1; 0 1;1 0; -1 0]

    mpQP = (H=H,f=f,F=F,
            A=A,b=b,B=B)
    Θ = (ub=1.5*ones(nth),lb=-1.5*ones(nth)) 

    opts = Dict("store_points"=>true, "verbose"=>1)
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    @test length(sol.CRs) == 27
end


@testset "Lower-dimensional critical regions with SCS" begin 
    n, nth = 2,1
    H = diagm(ones(n)) 
    f = zeros(n,1)
    F = zeros(n,nth)
    A = [0 -1.0;
         -1.0 0;
         0 -1;
         1 1;
         -1 0]
    b = [0;0;-1.0;3;-1]
    B = [[-0.5; -0.5; 0; 0; 0] zeros(5,0)]

    mpQP = (H=H,f=f,F=F,
            A=A,b=b,B=B)
    Θ = (ub=5*ones(nth),lb=1*ones(nth)) 

    opts = ParametricDAQP.Settings()
    opts.verbose=1;
    opts.store_points=true
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
end

@testset "Basic SISO Example from Bemporad et al. 2002" begin
    # Setup mpQP
    H =  [1.5064 0.4838; 0.4838 1.5258];
    f = zeros(2,1)
    F = [9.6652 5.2115; 7.0732 -7.0879];
    A = [1.0 0; -1 0; 0 1; 0 -1];
    B = zeros(4,2);
    b = 2*ones(4);
    mpQP = (H=H,f=f,F=F,A=A,b=b,B=B)

    # Setup parameter region of interest
    ub,lb  = 1.5*ones(2), -1.5*ones(2)
    Θ = (ub=ub,lb=lb)

    # Solve mpQP over desired region
    opts = ParametricDAQP.Settings()
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    bst = ParametricDAQP.build_tree(sol)

    # Compare bst with QP solution
    N = 1000
    ths = 1.5*(2*rand(2,N).-1);
    errs = zeros(N)
    for n = 1:N
        θ = ths[:,n]
        xbst = ParametricDAQP.evaluate(bst,θ./ub)
        f = mpQP.f[:,1]+mpQP.F*θ
        b = mpQP.b[:,1]+mpQP.B*θ
        xref,~,~,info= DAQP.quadprog(mpQP.H,f,mpQP.A,b);
        errs[n] = norm(xref-xbst)
    end
    @test all(errs.< 1e-5)
end

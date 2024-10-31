using ParametricDAQP
using PolyDAQP
using Test
using LinearAlgebra
using DAQP

function generate_mpQP(n,m,nth)
    M = randn(n,n)
    H = M*M'
    f = randn(n,1)
    f_theta = randn(n,nth)
    A = randn(m,n)
    b = [rand(m,1);rand(m,1)]
    F0 = randn(n,nth); # The point F0*th will be primal feasible
    W =[A;-A]*(-F0);
    #bounds_table = [collect(m+1:2m);collect(1:m)]
    bounds_table = collect(1:2m)
    senses = zeros(Cint,2m)
    mpQP = (H=H,f=f,f_theta=f_theta,H_theta=zeros(0,0),
                A=[A;-A],b=b,W=W,bounds_table=bounds_table,senses=senses)

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


    mpQP,Θ = generate_mpQP(n,m,nth)
    F,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);

    # Test for 10000 random points 
    N = 10000
    ths = 2*rand(nth,N).-1;
    containment_inds = zeros(Int,N) 
    errs = zeros(N) 
    for n = 1:N
        θ = ths[:,n]
        inds = pointlocation(θ,F;eps_gap=opts.eps_gap);
        containment_inds[n]=length(inds)
        xsol = F[inds[1]].xTH'*θ+F[inds[1]].xC
        f = mpQP.f[:,1]+mpQP.f_theta*θ
        b = mpQP.b[:,1]+mpQP.W*θ
        xref,~,~,info= DAQP.quadprog(mpQP.H,f,mpQP.A,b,-1e30*ones(2m),mpQP.senses);
        AS0 = findall(info.λ.!=0)
        errs[n] = norm(xsol-xref) 
    end
    @test all(errs.< tol)
    @test ~any(containment_inds.==0) # No holes
    @test sum(containment_inds.>1) < 0.01*N # Not more than 1% overlap
end

# Example 2 in Ahmadi-Moshkenani 2018
@testset "Full-dim LICQ" begin 
    n, nth = 4,2
    H = diagm(ones(n)) 
    f = zeros(n,1)
    f_theta = zeros(n,nth)
    A = [[1 0;
         -1 0;
         0 1;
         0 -1;] -ones(4,1) 0.5*ones(4,1)]
    b = -ones(4,1)
    W = [1.0 0; -1 0; 0 -1; 1 1];
    #bounds_table = [collect(m+1:2m);collect(1:m)]
    bounds_table = collect(1:length(b))
    senses = zeros(Cint,length(b))
    mpQP = (H=H,f=f,f_theta=f_theta,H_theta=zeros(0,0),
                A=A,b=b,W=W,bounds_table=bounds_table,senses=senses)
    Θ = (A = zeros(nth,0), b=zeros(0), ub=2*ones(nth),lb=-2*ones(nth)) 

    opts = ParametricDAQP.Settings()
    opts.verbose=1;
    opts.eps_gap = 0
    opts.postcheck_rank=true
    F,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
end

# Non facet to facet
@testset "Violated Facet-to-facet" begin 
    n, nth = 3,2
    H = diagm(ones(n)) 
    f = zeros(n,1)
    f_theta = zeros(n,nth)
    A = [1 0 -1;
         -1 0 -1;
         0 1 -1;
         0 -1 -1;
         3/4 16/25 -1;
         -3/4 -16/25 -1]
    b = -ones(6);
    W = [1.0 0; -1 0; 0 -1; 0 1;1 0; -1 0]
    #bounds_table = [collect(m+1:2m);collect(1:m)]
    bounds_table = collect(1:6)
    senses = zeros(Cint,6)
    mpQP = (H=H,f=f,f_theta=f_theta,H_theta=zeros(0,0),
                A=A,b=b,W=W,bounds_table=bounds_table,senses=senses)
    Θ = (A = zeros(nth,0), b=zeros(0), ub=1.5*ones(nth),lb=-1.5*ones(nth)) 

    opts = ParametricDAQP.Settings()
    opts.verbose=1;
    opts.store_points=true
    F,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    @test length(F) == 27
end


@testset "Lower-dimensional critical regions with SCS" begin 
    n, nth = 2,1
    H = diagm(ones(n)) 
    f = zeros(n,1)
    f_theta = zeros(n,nth)
    A = [0 -1.0;
         -1.0 0;
         0 -1;
         1 1;
         -1 0]

    b = [0;0;-1.0;3;-1]
    W = [[-0.5; -0.5; 0; 0; 0] zeros(5,0)]
    #bounds_table = [collect(m+1:2m);collect(1:m)]
    bounds_table = collect(1:5)
    senses = zeros(Cint,5)
    mpQP = (H=H,f=f,f_theta=f_theta,H_theta=zeros(0,0),
                A=A,b=b,W=W,bounds_table=bounds_table,senses=senses)
    Θ = (A = zeros(nth,0), b=zeros(0), ub=5*ones(nth),lb=1*ones(nth)) 

    opts = ParametricDAQP.Settings()
    opts.verbose=1;
    opts.store_points=true
    opts.eps_gap = 1e-7
    F,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
end

@testset "Basic SISO Example from Bemporad et al. 2002" begin
    # Setup mpQP
    H =  [1.5064 0.4838; 0.4838 1.5258];
    f = zeros(2,1)
    f_theta = [9.6652 5.2115; 7.0732 -7.0879];
    A = [1.0 0; -1 0; 0 1; 0 -1];
    W = zeros(4,2);
    b = 2*ones(4);
    mpQP = (H=H,f=f,f_theta=f_theta,A=A,b=b,W=W)

    # Setup parameter region of interest
    ub,lb  = 1.5*ones(2), -1.5*ones(2)
    Θ = (ub=ub,lb=lb)

    # Solve mpQP over desired region
    opts = ParametricDAQP.Settings()
    F,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
end

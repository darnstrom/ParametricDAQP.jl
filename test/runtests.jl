using ParametricDAQP
using PolyDAQP
using Test
using LinearAlgebra
using DAQPBase
const DAQP = DAQPBase
global templib

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
    rel_tol = 1e-5

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
        xsol = sol.CRs[inds[1]].z'*[θ;1]
        λsol = sol.CRs[inds[1]].lam'*[θ;1]
        f = mpQP.f[:,1]+mpQP.F*θ
        b = mpQP.b[:,1]+mpQP.B*θ
        xref,~,~,info= DAQP.quadprog(mpQP.H,f,mpQP.A,b,-1e30*ones(2m),mpQP.senses);
        @test xref ≈ xsol atol=1e-5 rtol=1e-5
        @test λsol ≈ info.λ[sol.CRs[inds[1]].AS] atol=1e-5 rtol=1e-5
    end
    @test ~any(containment_inds.==0) # No holes
    @test sum(containment_inds.>1) < 0.01*N # Not more than 1% overlap
end

@testset "Unsymmetric parameter region " begin
    n,m,nth = 3,10,4
    rel_tol = 1e-5

    opts = ParametricDAQP.Settings()
    opts.store_points=true
    opts.store_dual=true


    mpQP,Θ = generate_mpQP(n,m,nth)
    Θ.lb[:] = -(rand(nth).+2) # Make unsymmetric
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);

    # Test for 10000 random points 
    N = 10000
    ths = 2*rand(nth,N).-1;
    containment_inds = zeros(Int,N) 
    errs_primal,errs_dual = zeros(N),zeros(N)
    for n = 1:N
        θ = ths[:,n]
        CRs = ParametricDAQP.get_critical_regions(sol)
        inds = pointlocation(θ,CRs);
        containment_inds[n]=length(inds)
        AS = CRs[inds[1]].AS
        xsol = CRs[inds[1]].z'*[θ;1]
        λsol = CRs[inds[1]].lam'*[θ;1]
        f = mpQP.f[:,1]+mpQP.F*θ
        b = mpQP.b[:,1]+mpQP.B*θ
        xref,fval,exitflag,info= DAQP.quadprog(mpQP.H,f,mpQP.A,b,-1e30*ones(2m),mpQP.senses);
        @test xref ≈ xsol atol=1e-5 rtol=1e-5
        @test λsol ≈ info.λ[sol.CRs[inds[1]].AS] atol=1e-5 rtol=1e-5
    end
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

    opts = Dict("store_points"=>true, "verbose"=>1, "lowdim_tol" => 0)
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    @test length(sol.CRs) == 27
    # Remove lower dimensional regions
    opts = Dict("store_points"=>true, "verbose"=>1, "lowdim_tol"=>1e-9)
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    @test length(sol.CRs) == 12
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
        xbst = ParametricDAQP.evaluate(bst,θ)
        f = mpQP.f[:,1]+mpQP.F*θ
        b = mpQP.b[:,1]+mpQP.B*θ
        xref,~,~,info= DAQP.quadprog(mpQP.H,f,mpQP.A,b);
        errs[n] = norm(xref-xbst)
    end
    @test all(errs.< 1e-5)
end

@testset "Zero rows in A" begin
    n,m,nth = 3,10,4
    tol = 1e-5

    opts = ParametricDAQP.Settings()


    mpQP,Θ = generate_mpQP(n,m,nth)
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    Nref = length(sol.CRs)

    # Append some zero rows
    H = mpQP.H
    f = mpQP.f
    F = mpQP.F
    A = [mpQP.A;zeros(nth,n)] 
    b = [mpQP.b;ones(nth)]
    B = [mpQP.B;-I(nth)]
    bounds_table = [mpQP.bounds_table;2m+1:2m+nth]  
    senses = [mpQP.senses;zeros(Cint,nth)]  
    mpQP = (H=H,f=f,F=F,
                A=A,b=b,B=B,bounds_table=bounds_table,senses=senses)
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    @test Nref == length(sol.CRs)
end

@testset "Fixed parameters" begin
    n,m,nth = 3,10,4
    tol = 1e-5

    opts = ParametricDAQP.Settings()


    mpQP,Θ = generate_mpQP(n,m,nth)
    Θfix = (lb=[0;0;Θ.lb[3:4]],ub=[0;0;Θ.ub[3:4]])
    sol,info = ParametricDAQP.mpsolve(mpQP,Θfix;opts);
    Nref = length(sol.CRs)

    # Append some zero rows
    H = mpQP.H
    f = mpQP.f
    F = mpQP.F[:,3:4]
    A = mpQP.A
    b = mpQP.b
    B = mpQP.B[:,3:4]
    bounds_table = mpQP.bounds_table  
    senses = mpQP.senses
    mpQP = (H=H,f=f,F=F,
                A=A,b=b,B=B,bounds_table=bounds_table,senses=senses)
    sol,info = ParametricDAQP.mpsolve(mpQP,(lb=Θ.lb[3:4], ub = Θ.ub[3:4]);opts);
    @test Nref == length(sol.CRs)
end

@testset "Equality constraint" begin
    n,m,nth = 3,10,4
    tol = 1e-5

    mpQP,Θ = generate_mpQP(n,m,nth)
    mpQP = merge(mpQP,(eq_ids=[5],))

    opts = ParametricDAQP.Settings()
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    for cr in sol.CRs
        @test 5 ∈ cr.AS
    end
end

@testset "Primal degenerate LP" begin
    H =  zeros(4,4);
    f = [1.0;1.0;0;0] 
    F = zeros(4,2) 
    A = [-1.0 0 0 0;
         -1 0 -1 0;
         -1 0 0 0;
         -1 0 1 0;
         0 -1 -1 0;
         0 -1 -1 -1;
         0 -1 1 0;
         0 -1 1 1;
         0 0 1 0;
         0 0 -1 0;
         0 0 0 1;
         0 0 0 -1]
    B = [1.0 1;
         0 1;
         -1 -1;
         0 -1;
         1 2;
         0 1;
         -1 -2;
         0 -1;
         0 0;
         0 0;
         0 0;
         0 0]

    b = [zeros(8);ones(4)];
    mpQP = (H=H,f=f,F=F,A=A,b=b,B=B)

    # Setup parameter region of interest
    ub,lb  = 2.5*ones(2), -2.5*ones(2)
    Θ = (ub=ub,lb=lb)

    # Solve mpQP over desired region
    opts = ParametricDAQP.Settings()
    opts.lowdim_tol = 1e-9
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
end

@testset "Dual degenerate LP" begin
    H =  zeros(2,2);
    f = [-2.0;-1-0] 
    F = zeros(2,2) 
    A = [1.0 3;
         2 1;
         1 0;
         -1 0;
         0 -1]
    B = [-2.0 1;
         1 -2;
         1 1;
         0 0;
         0 0;
        ] 
    b = [9.0;8;4;0;0];

    mpQP = (H=H,f=f,F=F,A=A,b=b,B=B, out_inds=[1])

    # Setup parameter region of interest
    ub,lb  = 10*ones(2), -10*ones(2)
    Θ = (ub=ub,lb=lb)

    opts = ParametricDAQP.Settings()
    opts.lowdim_tol = 1e-9
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
end

# Test settings
@testset "Settings" begin
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
    opts.region_limit = 1;
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    @test sol.status == :RegionLimitReached
    opts.region_limit = 1e12;
    opts.time_limit = 0;
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    @test sol.status == :TimeLimitReached
    opts.region_limit = 1e12;
    opts.time_limit = 1e12;
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    @test sol.status == :Solved

    status = ParametricDAQP.codegen(sol; dir=tempname())
    @test status > 0
    status = ParametricDAQP.codegen(sol; dir=tempname(),max_reals = 1)
    @test status < 0
end


# Test settings
@testset "C-generated BST" begin
    # Setup mpQP
    H =  [1.5064 0.4838; 0.4838 1.5258];
    f = zeros(2,1)
    F = [9.6652 5.2115; 7.0732 -7.0879];
    A = [1.0 0; -1 0; 0 1; 0 -1];
    B = zeros(4,2);
    b = 2*ones(4);
    mpQP = (H=H,f=f,F=F,A=A,b=b,B=B)

    # Get a reference point
    θ = [-0.75;0.5]
    f = mpQP.f[:,1]+mpQP.F*θ
    b = mpQP.b[:,1]+mpQP.B*θ
    zref,~,~,info= DAQP.quadprog(mpQP.H,f,mpQP.A,b,-1e30*ones(4),zeros(Cint,4));
    λref = info.λ

    # Setup parameter region of interest
    ub,lb  = 1.5*ones(2), -1.5*ones(2)
    Θ = (ub=ub,lb=lb)

    opts = ParametricDAQP.Settings()
    opts.store_dual=true
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);

    srcdir = tempname()
    status = ParametricDAQP.codegen(sol; dir=srcdir)
    if(!isnothing(Sys.which("gcc")))
        testlib = "tree_test."* Base.Libc.Libdl.dlext
        run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib pdaqp.c`; dir=srcdir))
        z = zeros(Cfloat,2)
        global templib = joinpath(srcdir,testlib)
        ccall(("pdaqp_evaluate", templib), Cvoid, (Ptr{Cfloat}, Ptr{Cfloat}), Cfloat.(θ),z)
        @test norm(z - zref)/norm(z) < 1e-6
    end
    srcdir = tempname()
    status = ParametricDAQP.codegen(sol; dir=srcdir,dual=true)
    if(!isnothing(Sys.which("gcc")))
        testlib = "tree_test."* Base.Libc.Libdl.dlext
        run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib "pdaqp.c"`; dir=srcdir))
        z, λ = zeros(Cfloat, 2), zeros(Cfloat, 4)
        global templib = joinpath(srcdir,testlib)
        ccall(("pdaqp_evaluate", templib), Cvoid, (Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}), Cfloat.(θ),z,λ)
        @test norm(z - zref)/norm(z) < 1e-6
        @test norm(λ - λref)/norm(λ) < 1e-6
    end
end

@testset "Unconstrained" begin
    n,nth = 10,5
    mpQP,Θ = generate_mpQP(n,0,nth)
    opts = ParametricDAQP.Settings()
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);
    th = randn(nth)
    @test length(sol.CRs) == 1
    @test norm(sol.CRs[1].z'*[th;1] - (-mpQP.H\(mpQP.F*th +mpQP.f))) < 1e-6
end

@testset "Clipping for C-generated BST" begin
    # Setup mpQP
    H =  [1.5064 0.4838; 0.4838 1.5258];
    f = zeros(2,1)
    F = [9.6652 5.2115; 7.0732 -7.0879];
    A = [1.0 0; -1 0; 0 1; 0 -1];
    B = zeros(4,2);
    b = [1.0;2;3;4]
    mpQP = (H=H,f=f,F=F,A=A,b=b,B=B)

    # Setup parameter region of interest
    ub,lb  = 1.5*ones(2), -1.5*ones(2)
    Θ = (ub=ub,lb=lb)

    opts = ParametricDAQP.Settings()
    opts.store_dual=true
    sol,info = ParametricDAQP.mpsolve(mpQP,Θ;opts);

    srcdir = tempname()
    bst = ParametricDAQP.build_tree(sol)
    status = ParametricDAQP.codegen(sol; dir=srcdir, clipping=true)
    if(!isnothing(Sys.which("gcc")))
        testlib = "tree_test."* Base.Libc.Libdl.dlext
        run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib pdaqp.c`; dir=srcdir))
        z = zeros(Cfloat,2)
        global templib = joinpath(srcdir,testlib)
        for i = 1:25
            θ = 3*(rand(2).-0.5)
            f = mpQP.f[:,1]+mpQP.F*θ
            b = mpQP.b[:,1]+mpQP.B*θ
            zref,~,~,info= DAQP.quadprog(mpQP.H,f,mpQP.A,b,-1e30*ones(4),zeros(Cint,4));
            ccall(("pdaqp_evaluate", templib), Cvoid, (Ptr{Cfloat}, Ptr{Cfloat}), Cfloat.(θ),z)
            @test norm(z - zref)/norm(z) < 1e-6
            @test norm(ParametricDAQP.evaluate(bst,θ) - zref)/norm(zref) < 1e-6
        end
    end
end

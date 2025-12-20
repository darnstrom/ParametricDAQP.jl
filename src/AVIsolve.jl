@enum WarmStartType NoWarmStart UnconstrainedSolution

function AVIsolve(
    M::AbstractMatrix,
    q::AbstractVector,
    A::AbstractMatrix,
    b::AbstractVector;
    max_iter=10^5,
    stepsize=0.5,
    tol=10^(-5),
    warmstart::WarmStartType=NoWarmStart)

    # Solves VI(Hx +h, C) where C={x| Ax <= b}
    # Using Douglas-Rachford
    m,n = size(A)
    M1 = (M + M') / 4
    M2 = M1 + (M - M') / 2
    M1plusH = M1 + I
    M2minusH = M2 - I
    M2plusH = lu!(M2 + I)
    if warmstart == NoWarmStart
        x = zeros(n)
    elseif warmstart == UnconstrainedSolution
        x = -M \ q
    end
    # check for trivial infeasibility: A[i,:] = 0, b[i]<0 for some i
    zero_rows = findall(row -> norm(row) < tol, eachrow(A))
    for idx_row in zero_rows
        b[idx_row] < 0 && return nothing, NaN, :Infeasible
    end
    # MAIN LOOP
    M2x_plus_q, y = zeros(n),zeros(n)
    backstep = DAQPBase.Model()
    DAQPBase.setup(backstep, M1plusH, M2x_plus_q, A, b, Float64[], zeros(Cint, m))
    proj = DAQPBase.Model()
    DAQPBase.setup(proj, Matrix{Float64}(I, n, n), -y, A, b,Float64[], zeros(Cint, m))
    for k in 1:max_iter
        M2x_plus_q[:] = q
        mul!(M2x_plus_q, M2minusH, x, 1.0, 1.0)
        DAQPBase.update(backstep, nothing, M2x_plus_q, nothing, nothing,nothing)
        y[:], _, exitflag, _ = DAQPBase.solve(backstep)
        if exitflag < 0
            return nothing, NaN, :Infeasible # infeasible
        end
        if isempty(y)
            error("[AVIsolve] Error when performing the projection step at iteration $k")
        end
        y[:] = 2 * stepsize * y + (1 - 2 * stepsize) * x
        mul!(y, M2, x, 1.0, 1.0)
        ldiv!(x, M2plusH, y)

        # Compute residual
        if mod(k, 10) == 0
            y = x - (M * x + q)
            DAQPBase.update(proj, nothing, -y, nothing, nothing,nothing)
            x_transf, _, _, _ = DAQPBase.solve(proj)
            r = norm(x - x_transf)
            # println("[AVIsolve]: Iteration = $k; Residual = $r")
            if r < tol
                return x, r, :Solved
            elseif k == max_iter
                @warn "[AVIsolve.jl] Maximum iterations reached, residual = $r"
                return x, r, :MaximumIterationsReached
            end
        end
    end

end

# Fixed point iteration for implicit methods
function FixIter(fun, u; tol=sqrt(eps(Float64)), max_iter=1000)
    for _ in 1:max_iter
        uold = u
        u = fun(u)
        norm(uold - u) ≤ tol && return u
    end
    error("FixIter failed to converge after $(max_iter) iterations (tol=$(tol)).")
end

# First order explicit Euler (not symplectic)
function ExpEuler(f, z, h)
    return z + h * f(z)
end

# First order implicit Euler (not symplectic)
function ImpEuler(f, z, h)
    g(k) = f(z + h * k)
    k = FixIter(g, f(z); tol=1e-12)
    return z + h * k
end

# Symplectic Euler for separable Hamiltonians: f1(q,p) = dq/dt, f2(q,p) = dp/dt
function SympEuler(f1, f2, z, h)
    n = Int(length(z) / 2)
    q = z[1:n]
    p = z[n+1:end]

    p_new = p + h * f2(q, p)
    q_new = q + h * f1(q, p_new)

    return [q_new; p_new]
end

# 4th order Runge-Kutta method (not symplectic)
function RK4(f, z, h)
    K1 = f(z)
    K2 = f(z + h * K1 / 2)
    K3 = f(z + h * K2 / 2)
    K4 = f(z + h * K3)
    return z + h * (K1 + 2 * K2 + 2 * K3 + K4) / 6
end

# Implicit midpoint rule (symplectic)
function MidPoint(f, z, h)
    g(K) = f(z + h * K / 2)
    K0 = f(z)
    K = FixIter(g, K0; tol=1e-12)
    return z + h * K
end

# Stormer-Verlet method for separable Hamiltonians: f1(q,p) = dq/dt, f2(q,p) = dp/dt
function StormerVerlet(f1, f2, z, h)
    n = Int(length(z) / 2)
    q0 = z[1:n]
    p0 = z[n+1:end]

    g1(Km) = q0 + (h / 2) * f1(Km, p0)
    Km0 = q0 + (h / 2) * f1(q0, p0)
    Kmid = FixIter(g1, Km0; tol=1e-14)

    g2(K) = p0 + (h / 2) * (f2(Kmid, p0) + f2(Kmid, K))
    K0 = p0
    Knew = FixIter(g2, K0; tol=1e-14)

    p = Knew
    q = Kmid + (h / 2) * f1(Kmid, Knew)
    return [q; p]
end

# Generic integrator: applies numFlow repeatedly for the given number of steps
function Integrator(numFlow, z0, steps)
    d = length(z0)
    z = zeros(d, steps + 1)
    z[:, 1] = z0
    for k = 1:steps
        z[:, k+1] = numFlow(z[:, k])
    end
    return z
end

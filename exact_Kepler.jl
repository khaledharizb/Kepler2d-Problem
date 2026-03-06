# From textbook (auth: Bruno Cordani): The Kepler Problem pages:25-27

function phiKepler(z, h, ecc::Float64)
    tol = 1e-12
    q₀ = z[1:2]
    p₀ = z[3:4]
    r₀ = norm(q₀)
    H₀ = 0.5 * p₀' * p₀ - 1 / r₀
    a = -1 / (2 * H₀)
    w = a^(3 / 2)

    g(E) = ecc * sin(E) + h / w                            # eq. (2.2.6)
    E = FixIter(g, w * h; tol=tol)

    q = [a * cos(E) - a * ecc;
         a * sqrt(1 - ecc^2) * sin(E)]                     # page 27

    dtdE = w * (1 - ecc * cos(E))                           # line before eq. (2.2.6)
    p = [-a * sin(E) / dtdE;
          a * sqrt(1 - ecc^2) * cos(E) / dtdE]             # p = dq/dt

    return [q; p]
end

function phiKepler(z, h, steps::Int, ecc::Float64)
    q = zeros(length(z), steps + 1)
    q[:, 1] = z
    t0 = 0.0
    for i = 1:steps
        t0 += h
        q[:, i+1] = phiKepler(z, t0, ecc)
    end
    return q
end

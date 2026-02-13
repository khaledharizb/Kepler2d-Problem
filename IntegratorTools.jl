# Fixed point iteration for implicit methods
function FixIter(fun, u; tol=sqrt(eps(Float64)), max_iter=1000) 
    for _ in 1:max_iter
        uold = u
        u = fun(u)
        norm(uold - u) ≤ tol && return u
    end
error("FixIter failed to converge after $(max_iter) iterations (tol=$(tol)).")
end

function  ExpEuler(f,z,h) # first order explicit Euler (not symplectic)
   return z + h * f(z); 
end;

function  ImpEuler(f,z,h) # first order implicit Euler (not symplectic)
 g(k) =  f(z + h * k);
 k = FixIter(g,f(z),1e-12);
 return z + h * k;   
end;

function SymEuler(f, z, h)
    d = length(z);
    m = Int(d / 2);

    znew = copy(z)

    F = f(z);

    # update q in-place in znew
    @inbounds for i = 1:m
        znew[i] += h * F[i]
    end

    # compute new vector field at mixed state
    Ftemp = f(znew)

    # update p
    @inbounds for i = m+1:d
        znew[i] += h * Ftemp[i]
    end

    return znew
end



function  RK4(f,z,h)  # 4th order RK method (not symplectic)
  K1 = f(z)     
  K2 = f(z + h * K1 / 2);
  K3 = f(z + h * K2 / 2);
  K4 = f(z + h * K3);
 return  z + h * (K1 + 2 * K2 + 2 *  K3 + K4) / 6;  
end;

function  MidPoint(f,z,h) # Midpoint rule (symplectic)
 g(K) =   f( z + h * K / 2 );
   K0 = f(z)
    K = FixIter(g,K0,1e-12);
return  z + h * K;
end;
  
function  StormerVerlet(f1,f2,z,h) # Stormer Verlet methods fur non-separable Hamiltonian f1=f(q,p), f2=g(q,p)

n = Int(length(z)/2); 
q0 = z[1:n];
p0 = z[n+1:end];

    
g1(Km) = q0 + (h/2) * f1(Km,p0);
      Km0 = q0 + (h/2) * f1(q0,p0); 
      Kmid = FixIter(g1,Km0,1e-14);
        
g2(K) = p0 + (h/2) * (f2(Kmid,p0)+f2(Kmid,K));
        K0 = p0; 
      Knew = FixIter(g2,K0,1e-14);

   p = Knew;    
   q = Kmid + (h/2) * f1(Kmid,Knew);
  return  return [q; p]; 
end;

function Integrator(numFlow,z0,steps)
    d = length(z0);
    z = zeros(d,steps+1)
    z[:,1] = z0
    for k = 1:steps
        z[:,k + 1] = numFlow(z[:,k])
    end
    return z
end;

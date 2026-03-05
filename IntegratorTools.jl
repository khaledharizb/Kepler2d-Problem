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
 k = FixIter(g, f(z); tol=1e-12);
 return z + h * k;   
end;

function SympEuler(f1,f2,z,h) # Symplectic Euler for non-separable Hamiltonian f1=f(q,p), f2=g(q,p)

 n = Int(length(z)/2);
q0 = z[1:n];
p0 = z[n+1:end];
        
g(K) = p0 + h * f2(q0,K);
  K0 = p0; 
   K = FixIter(g, K0; tol=1e-12);       
   p = K;   
   q = q0 + h * f1(q0,p);               
  return [q; p]; 
end;

function SympEuler(f,z,h) # Symplectic Euler using a single vector field f(q,p) -> [qdot; pdot]
 n = Int(length(z)/2)
 q0 = z[1:n]
 p0 = z[n+1:end]

 function qdot(q,p)
    return f(q,p)[1:n]
 end

 function pdot(q,p)
    return f(q,p)[n+1:end]
 end

 g(K) = p0 + h * pdot(q0, K)
 p = FixIter(g, p0; tol=1e-12)
 q = q0 + h * qdot(q0, p)
 return [q; p]
end;


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
    K = FixIter(g, K0; tol=1e-12);
return  z + h * K;
end;
  
function  StormerVerlet(f1,f2,z,h) # Stormer Verlet method for non-separable Hamiltonian f1=f(q,p), f2=g(q,p)

n = Int(length(z)/2); 
q0 = z[1:n];
p0 = z[n+1:end];

    
g1(Km) = q0 + (h/2) * f1(Km,p0);
      Km0 = q0 + (h/2) * f1(q0,p0); 
      Kmid = FixIter(g1, Km0; tol=1e-14);
        
g2(K) = p0 + (h/2) * (f2(Kmid,p0)+f2(Kmid,K));
        K0 = p0; 
      Knew = FixIter(g2, K0; tol=1e-14);

   p = Knew;    
   q = Kmid + (h/2) * f1(Kmid,Knew);
  return [q; p]; 
end;

function Integrator(numFlow,z,steps)
    Z = zeros(4,steps+1)
    Z[:,1] = z
    for k = 1:size(Z,2)-1
        Z[:,k + 1] = numFlow(Z[:,k])
    end
    return Z
end;

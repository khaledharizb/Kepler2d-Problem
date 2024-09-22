function FixIter(fun,u,tol) # Fixed point iteration used for implicit methods
    while true
        uold = u;
        u = fun(u);
        norm(uold - u) > tol || break;
    end
    return u;
end

function  ExpEuler(f,z,h) # first order explicit Euler (not symplectic)
   return z + h * f(z); 
end;

function  ImpEuler(f,z,h) # first order implicit Euler (not symplectic)
 g(k) =  f(z + h * k);
 k = FixIter(g,f(z),1e-12);
 return z + h * k;   
end;

function SympEuler(f1,f2,z,h) # Symplectic Euler fur non-separable Hamiltonian f1=f(q,p), f2=g(q,p)

 n = Int(length(z)/2);
q0 = z[1:n];
p0 = z[n+1:end];
        
g(K) = p0 + h * f2(q0,K);
  K0 = p0; 
   K = FixIter(g,K0,1e-12);       
   p = K;   
   q = q0 + h * f1(q0,p);               
  return [q; p]; 
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

function Integrator(numFlow,z,steps)
    Z = zeros(4,steps+1)
    Z[:,1] = z
    for k = 1:size(Z,2)-1
        Z[:,k + 1] = numFlow(Z[:,k])
    end
    return Z
end;
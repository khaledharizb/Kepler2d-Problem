# From textbook (auth: Bruno Cordani): The Kepler Problem pages:25-27 

function phiKepler(z,h) 
tol=1.e-12; q₀ = z[1:2] ; p₀ = z[3:4] ; r₀ = norm(q₀);     H₀ = 0.5 * p₀' * p₀ - 1 / r₀ ; a = -1 / (2*H₀) ; w = a^(3/2) ;
g(E) = ECC * sin(E) + h / w ;
  x₀ = w * h;  
   E = FixIter(g,x₀,tol);
   q = [a * cos(E) - a * ECC ; a * sqrt(1-ECC^2) * sin(E)]; 
dsdt = w * (1 - ECC * cos(E));    
   p = [- sin(E) / dsdt; a * sqrt(1-ECC^2) * cos(E) / dsdt];
  return [q;p];
end;

function phiKepler(z,h,steps) 
q = zeros((size(z)[1],steps+1));
q[:,1] = z;    
t0  = 0;
for i=1:steps;
  t0 = t0+h;
 q[:,i+1] = phiKepler(z,t0) ; 
end   
    return q;
end;
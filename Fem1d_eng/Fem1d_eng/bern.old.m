function [bp,bn]=bern(x)
%
% Funzione per il calcolo accurato della Bernoulli function
%
xlim=1e-2; 
ax=abs(x);
if (ax == 0)
   bp=1.; 
   bn=1.; 
return
end
if (ax > 80),
   if (x > 0)
     bp=0.; 
     bn=x;  
     return
   else      
     bp=-x; 
     bn=0.; 
     return  
   end
end
if (ax > xlim)
   bp=x/(exp(x)-1); 
   bn=x+bp; 
   return
else
   ii=1; 
   fp=1.;
   fn=1.; 
   df=1.; 
   s=1.;
   while (abs(df) > eps),
     ii=ii+1;  
     s=-s; 
     df=df*x/ii;
     fp=fp+df; 
     fn=fn+s*df;
     bp=1./fp; 
     bn=1./fn;
   end
   return
end
return

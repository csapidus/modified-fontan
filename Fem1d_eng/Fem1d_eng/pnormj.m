function [fnorm] =pnormj(n,alpha,beta) 
% 
% 
% 
const = alpha+beta+1; 
if n <= 1  
  pa = gammaf(n+alpha); 
  pb = gammaf(n+beta); 
  pn = gammaf(n); 
  pnab = gammaf(n+alpha+beta); 
  fnorm = (pa*pb)/(pn*pnab)*2^const/(2*n+const); 
  return 
end 
pa = gammaf(1+alpha); 
pb = gammaf(1+beta); 
pc = gammaf(1+const); 
prod = (pa*pb)/(2*(1+const)*pc)*(1+alpha)*(2+alpha)*... 
        (1+beta)*(2+beta); 
for i =3:n 
  frac = (i+alpha)*(i+beta)/(i*(i+alpha+beta)); 
  prod = prod*frac; 
end 
fnorm = prod*2^const/(2*n+const); 
% 
return

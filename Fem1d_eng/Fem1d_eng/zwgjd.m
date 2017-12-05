function [x,w] = zwgjd(np,alpha,beta) 
%  
%ZWGJD  Gauss-Legendre-Lobatto quadrature formula.
%
%    [X,W]=XWGJD(NP,ALPHA,BETA) returns the interior nodes 
%    X in (-1,1) and the related weights W of the 
%    Gauss-Legendre-Lobatto quadrature formula
%    with NP nodes. 
%
%    See also, XWGLL.  
%    Requirements: JACG, JACOBF, PNORMJ.

%   Gervasio Paola, 1997. Saleri Fausto, 1999.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:48 $

n   = np-1; 
apb = alpha+beta; 
if np == 1 
 x(1) = (beta-alpha)/(apb+2); 
 g1 = gammaf(alpha+1); 
 g2 = gammaf(beta+1); 
 g3 = gammaf(apb+2); 
 w(1) = (g1*g2/g3)*2^(apb+1); 
 return 
end 
% 
x = jacg(np,alpha,beta); 
np1 = n+1; 
np2 = n+2; 
fac1 = np1+alpha+beta+1; 
fac2 = fac1+np1; 
fac3 = fac2+1; 
fnorm =pnormj(np1,alpha,beta); 
rcoef = (fnorm*fac2*fac3)/(2*fac1*np2); 
for i=1:np 
 [p,pd,pm1,pdm1,pm2,pdm2]= jacobf(np2,alpha,beta,x(i)); 
 w(i) = -rcoef/(p*pdm1); 
end 
 
return

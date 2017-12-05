function [nqn,qnodes,w] = gausscheb(np) 
%
%GAUSSCHEB    One-dimensional Gauss-Chebyshev quadrature rule.
%
%    [NQN,QNODES,W]=GAUSSCHEB(NP) returns the
%    NP quadrature nodes QNODES and the NP weights
%    of the Gauss-Chebysev quadrature rule in (0,1).
 
%   Gervasio Paola, 1997.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:47 $

nqn = np;
den=pi/(2*np); 
ww=pi/np; 
for j = 0:np-1    
   qnodes(j+1) = -cos((2*j+1)*den);    
   w(j+1) = ww; 
end 

qnodes = 1/2*(qnodes+1);
w = w/2;
   
return
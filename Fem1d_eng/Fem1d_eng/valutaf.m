function [nq,bq,sq,fq] = valutaf (nqn, qnodes, h, xi,...
                                  nu, beta, sigma, f,t)
%
%VALUTAF        Evaluation of the coefficient of an elliptic
%               problem in the quadrature nodes. 
%
%    [NQ,BQ,SQ,FQ]=VALUTAF(NQN,QNODES,H,XI,MU,BETA,SIGMA,F) 
%    return in the array NQ, BQ, SQ and FQ the evaluation of
%    the coefficients NU, BETA, SIGMA and F in the quadrature
%    nodes QNODES.
%
 
%   Saleri Fausto, 1999.

if nargin == 6 | nargin==7
sq = 0; fq = 0;if nargin ==7, t=sigma; end;
for k = 1:nqn
   x = h*qnodes(k)+xi;
   nq (k) = eval(nu);
   bq(k)  = eval(beta);
end
else
for k = 1:nqn   
   x = h*qnodes(k)+xi;
   fq (k) = eval(f);
   nq (k) = eval(nu);
   bq (k) = eval(beta);
   sq (k) = eval(sigma);   
end
end

return

function [nqn, qnodes, w] = quadratura (fem)
%
%   [nqn, qnodes, w] = quadratura (fem)
%
% Formula di quadratura  sull'intervallo (0,1) 
%
% Se fem=1, viene usata la formula di Simpson
% (grado di precisione 2), in caso contrario una
% formuladi quadratura di Gauss-Legendre con grado
% di precisione 5.
%
%if fem == 1
%   nqn = 3;
%   qnodes = [0 0.5 1];
%   w = [1/6 2/3 1/6];
%else
   nqn = 3;
   qnodes = [0.11270166537926 0.50000000000000 0.88729833462074];
   w = [0.55555555555556 0.88888888888889 0.55555555555556]*0.5;
%end
return

function [aloc,bloc] = locassembla_nc (nln, nqn, mu, sigma, beta, ...
   h, w, dphiq, phiq, fq)
%
%LOCASSEMBLA    Local matrix and local load vector for an
%               elliptic one-dimensional problem. 
%
%    [PHIQ,DPHIQ]=VALUTA(FEM,NQN,QNODES,PHI,DPHI) 
%    return in the array PHIQ and DPHIQ the evaluation of
%    the basis function and their derivatives in the quadrature
%    nodes QNODES.
%
 
%   Saleri Fausto, 1999.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:48 $


aloc = zeros (nln);
bloc = zeros (nln,1);

for i = 1:nln
   for j = 1:nln
      aloc(i,j) = sum(w'.*((1/h)*mu.*dphiq(i,:).*dphiq(j,:)+...         
      h*sigma.*phiq(i,:).*phiq(j,:)+...         
      beta.*phiq(i,:).*dphiq(j,:)));
   end   
end
for i = 1:nln
   bloc (i) = sum(h*w'.*fq.*phiq(i,:));
end

return

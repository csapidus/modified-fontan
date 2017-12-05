function [errore]=h1errore (fem, ne, ndof, nln, dof, list_el, ufem, u, ux,t)
%
%  Calcolo dell'errore in norma H1
%
[nqn, qnodes, w] = gausslegendre(7);
[phi,dphi] = bases1d (fem);
errL2  = 0;
errH1  = 0;
errsH1 = 0;
%figure(1); hold on;
for ie = 1:ne
   jac = ( dof (list_el(2,ie)) - dof (list_el(1,ie)) );
   [phiq, dphiq] = valuta (fem, nqn, qnodes, phi, dphi);
   if nargin==10, [uq, uxq] = valutaf (nqn, qnodes, jac, dof(list_el(1,ie)), u, ux,t);
   else, [uq, uxq] = valutaf (nqn, qnodes, jac, dof(list_el(1,ie)), u, ux);end;   
   uhq  = zeros(nqn,1);
   uhxq = zeros(nqn,1);
   for iloc = 1:nln
      i = list_el(iloc,ie);
      for k = 1:nqn
         uhxq (k) = uhxq(k)+ufem(i)*dphiq(iloc,k)/jac;
         uhq  (k) = uhq(k)+ufem(i)*phiq(iloc,k);
      end
   end
   x = jac*qnodes+dof(list_el(1,ie));
%   plot(x,uhxq)
   for k = 1:nqn
      errsH1 = errsH1 + w(k)*jac*(uxq(k)-uhxq(k))^2;
      errL2  = errL2  + w(k)*jac*(uq(k)-uhq(k))^2;
   end
end
%hold off
errore(1) = sqrt(errL2 + errsH1);
errore(2) = sqrt(errsH1);
errore(3) = sqrt(errL2);

return

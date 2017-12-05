function [M,u0]=crea_massa(fem,ne,ndof,nln,dof,list_el,u_init)
%
% Creazione della matrice di Massa (pb parabolico)
% e del termine noto per il proiettato del dato iniziale
% A. Veneziani - Aprile 2000
%
[nqn, qnodes, w] = xwgll(7)

[phi,dphi] = bases1d (fem);

M = sparse (ndof,ndof); 
u0 = zeros(ndof,1);
for ie = 1:ne   
   jac= dof (list_el(2,ie)) - dof (list_el(1,ie));  
   [phiq, dphiq] = valuta (fem, nqn, qnodes, phi, dphi);
   if (nargin>6), u0q = valutaf(nqn,qnodes,jac,dof(list_el(1,ie)),...
         u_init,u_init,u_init,u_init); end;  
   Mloc = zeros (nln);

   for i = 1:nln
    for j = 1:nln
      Mloc(i,j) = jac*sum(w'.*(phiq(i,:).*phiq(j,:)));
    end   
   end
   if (nargin>6)
    for i = 1:nln
       u_loc(i) = jac*sum(w'.*u0q.*phiq(i,:));
    end
   end
     
 %
   for iloc = 1:nln
      i = list_el (iloc,ie);
      for jloc = 1:nln
         j = list_el (jloc,ie);
         M (i,j) = M (i,j) + Mloc(iloc,jloc);
      end
   end
   if (nargin>6)
    for iloc = 1:nln
      i = list_el (iloc,ie);
      u0(i) = u0(i) + u_loc(iloc);
    end
   end  
end
return
   
      
   
   
   
   
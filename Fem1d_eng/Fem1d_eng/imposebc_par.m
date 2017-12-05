function [uh,matout,bout,stima_cond]=imposebc_par(matfem,bfem,ndof,bc,massa,matold,bold,...
                                      u_old,theta,DeltaT,types)
%
%   Veneziani Alessandro 2000
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:48 $

if nargin == 10
   types = 0;
end

if length(bc)==2
   %
   % Dirichlet boundary conditions: mettere i rilevamenti SOLO su t_noto
   bfem = bfem - matfem(:,1)*bc(1) - matfem(:,ndof)*bc(2);
   bout = bfem([2:ndof-1]);
   matout = matfem([2:ndof-1],[2:ndof-1]);
   aux = 1/DeltaT*(massa(:,1)*bc(1) + massa(:,ndof)*bc(2));
   aux2 = 1/DeltaT*massa*u_old;
   %
   matrice = 1/DeltaT*massa([2:ndof-1],[2:ndof-1])+theta*matout;
   t_noto= aux2(2:ndof-1) + ...
      theta*bout + (1-theta)*(-matold*u_old(2:ndof-1)+bold) - aux(2:ndof-1);
   %   
   [ufem]=fsolver(matrice,t_noto,types);
   uh = [bc(1); ufem; bc(2)];
elseif bc(3)>0 
   %   
   % Dirichlet type in the last node  <=> Robin/Neumann in the first
   bfem = bfem - matfem (:,ndof)*bc(2);
   bfem(1) = bfem(1) + bc(1);
   bout = bfem([1:ndof-1]);
   matfem(1,1) = matfem(1,1) + bc(4);
   matout = matfem([1:ndof-1],[1:ndof-1]);
   %
   matrice = 1/DeltaT*massa([1:ndof-1],[1:ndof-1])+theta*matout;
   aux=massa*u_old;
   t_noto=1/DeltaT*aux(1:ndof-1) + theta*bout + (1-theta)*(matold*u_old(1:ndof-1)+bold);
   %   
   [ufem]=fsolver(matrice,t_noto,types);
   uh = [ufem;bc(2)];
elseif bc(3)<0
   %
   % Dirichlet type in the first node <=> Robin/Neumann in the last
   bfem = bfem(ndof) + bc(2);
   bfem = bfem - matfem (:,1)*bc(1);
   bout = bfem([2:ndof]);
   matfem(ndof,ndof) = matfem(ndof,ndof) + bc(5);
   matout = matfem([2:ndof],[2:ndof]);
   %
   matrice = 1/DeltaT*massa([2:ndof],[2:ndof])+theta*matout;
   aux=massa*u_old;
   t_noto=1/DeltaT*aux(2:ndof) + theta*bout + (1-theta)*(matold*u_old(1:ndof-1)+bold);
   %   
   [ufem]=fsolver(matrice,t_noto,types);
   uh = [bc(1); ufem];
else
   %
   % Robin (hopefully!) on both the endpoints
   bfem(1) = bfem(1) + bc(1);
   bfem(ndof) = bfem(ndof) + bc(2);
   matfem(1,1) = matfem(1,1) + bc(4);
   matfem(ndof,ndof) = matfem(ndof,ndof) + bc(5);
   bout=bfem; matout=matfem;
   %
   matrice = 1/DeltaT*massa+theta*matout;
   t_noto=1/DeltaT*massa*u_old + theta*bout + (1-theta)*(matold*u_old(1:ndof)+bold);
   %   
   [uh]=fsolver(matfem,bfem,types);
end;
   stima_cond=condest(matrice);

return
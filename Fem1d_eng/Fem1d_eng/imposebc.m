function [uh,matfem]=imposebc(matfem,bfem,ndof,bc,types,mu_sx,mu_dx)
%
%IMPOSEBC       Assign the boundary conditions and solve
%               the finite element linear system.
%
%    [MATFEM,BFEM]=IMPOSEBC(MATFEM,BFEM,NDOF,BC) impose
%    in the matrix MATFEM and in the load vector BFEM the
%    boundary conditions. The function SOLVSYS solve the
%    associated linear system with a direct method.
%    If bc is a 1x2 array, the boundary conditions are
%    of Dirichlet type at the extrema of Omega:
%        u(Omega(1)) = bc(1),      u(Omega(2)) = bc(2);
%    if bc is a 1x3 array, the boundary conditions are:
%    bc(3) < 0 of Neumann type in Omega(1), of Dirichlet 
%    type in Omega(2) and
%        u'(Omega(1)) = bc(1),  u(Omega(2)) = bc(2);
%    bc(3) > 0 of Dirichlet type in Omega(1), of Neumann 
%    type in Omega(2) and
%        u(Omega(1)) = bc(1),  u'(Omega(2)) = bc(2);
%    bc(3) = 0 of Neumann type in Omega(1) and Omega(2) and
%        u'(Omega(1)) = bc(1),  u'(Omega(2)) = bc(2).
%    In this case sigma should be not equal zero.
%
%     [MATFEM,BFEM]=IMPOSEBC(MATFEM,BFEM,NDOF,BC,TYPES)
%    If TYPES is equal 0 (default) a direct method is used, 
%    otherwise an iterative method will be used (see, 
%
 
%   Saleri Fausto, 1999.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:48 $

if nargin == 4
   types = 0;
end

if length(bc)==2
   %
   % Dirichlet boundary conditions
   bfem = bfem - matfem (:,1)*bc(1) - matfem (:,ndof)*bc(2);
   bfem = bfem([2:ndof-1]);
   matfem = matfem([2:ndof-1],[2:ndof-1]);
   [ufem]=fsolver(matfem,bfem,types);
   uh = [bc(1); ufem; bc(2)];
elseif bc(3) < 0
   %   
   % Dirichlet type in the last node  <=> Robin/Neumann in the first
   bfem = bfem - matfem (:,ndof)*bc(2);
   bfem(1) = bfem(1) + bc(1)*mu_sx;
   bfem = bfem([1:ndof-1]);
   matfem(1,1) = matfem(1,1) + bc(4)*mu_sx;
   matfem = matfem([1:ndof-1],[1:ndof-1]);
   [ufem]=fsolver(matfem,bfem,types);
   uh = [ufem;bc(2)];
elseif bc(3) > 0
   %
   % Dirichlet type in the first node <=> Robin/Neumann in the last
   bfem(ndof) = bfem(ndof) + bc(2)*mu_dx;
   bfem = bfem - matfem (:,1)*bc(1);
   bfem = bfem([2:ndof]);
   matfem(ndof,ndof) = matfem(ndof,ndof) + bc(5)*mu_dx;
   matfem = matfem([2:ndof],[2:ndof]);
   [ufem]=fsolver(matfem,bfem,types);
   uh = [bc(1); ufem];
else
   %
   % Robin (hopefully!) on both the endpoints
   bfem(1) = bfem(1) + bc(1)*mu_sx;
   bfem(ndof) = bfem(ndof) + bc(2)*mu_dx;
   matfem(1,1) = matfem(1,1) + bc(4)*mu_sx;
   matfem(ndof,ndof) = matfem(ndof,ndof) + bc(5)*mu_dx;
   [uh]=fsolver(matfem,bfem,types);
end

return
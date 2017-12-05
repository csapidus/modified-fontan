function [dof,ndof,nln,flag] = create_dof (fem,coord,nvert,flag)
%   
%CREATE_DOF     Degree of freedom of one-dimensional 
%               finite elements. 
%
%    [DOF,NDOF,NLN]=CREATE_DOF(FEM,COORD,NVERT) return in the
%    1 x NDOF array DOF the coordinates of the degree of freedom
%    of the finite element FEM. NLN is the local number of degree
%    of freedom.
%
 
%   Saleri Fausto, 1999.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:47 $

% Check for an acceptable type of finite element

switch nargin
case 4,   
   if flag < 0, dof = []; ndof = []; return; end
otherwise,
   switch nargout
   case 4,
      flag = 0;
   end
end

switch abs(fem)
case 1,
   dof = coord;   ndof = nvert;  nln = 2;
case 2,
   k = 0; dof = [];
   for i = 1:nvert-1
      xm = (coord(i+1)+coord(i))*0.5;
      dof = [dof, coord(i),xm];
      k = k + 2;
   end
   ndof = k+1;  nln = 3;
   dof (ndof) = coord (nvert);
case 3,
   k = 0; dof = [];
   for i = 1:nvert-1
      x13 = (coord(i+1)+2*coord(i))/3;
      x23 = (2*coord(i+1)+coord(i))/3;
      dof = [dof, coord(i),x13,x23];
      k = k + 3;
   end
   ndof = k+1;  nln = 4;
   dof (ndof) = coord (nvert);
otherwise,
   es=sprintf(['Finite elements available are P1, P2 and P3']);
   disp(es); flag=-1; dof=[]; ndof=0; nln=0; return;
end

return
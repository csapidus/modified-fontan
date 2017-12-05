function [uh,coord,out1,out2]=fem1d_ell_c(Omega,bc,femsys,coeff,hmesh,arg1,arg2,arg3)
%   
%FEM1D        One-dimensional finite element solver for linear 
%             elliptic problem with non constant coefficients.
%
%    [UH,COORD] = FEM1D(OMEGA,BC,FEMSYS,COEFF,HMESH) returns the
%   finite element solution UH of the elliptic problem in the interval
%   (OMEGA(1),OMEGA(2)):
%        -(NU(x) u'(x))' + (BETA(x) u(x))' + SIGMA u(x) = F(x)  
%   For the format of the coefficients in the string COEFF see
%   the function DATAEXTR.
%    FEMSYS(1) is a variable equal to k for Pk elements, k=1,2,3.
%   Negative values of FEMSYS(1) denotes the corresponding yerarchical
%   elements.
%    FEMSYS(2) is a variable to set the solver for the 
%   finite element linear system (see, FSOLVER).
%    HMESH is a macro for the grid stepsize. For uniform grid
%   HMESH is a positive number. 
%   In the output variable COORD is stored the computed grid.
%
%    [UH,COORD] = FEM1D(OMEGA,BC,FEMSYS,COEFF,HMESH,ADAPT) realize
%   a simple grid adaption in order to reduce the H1-error tunder
%   the positive tolerance ADAPT.
%
%    [UH,COORD,ERROR] = FEM1D(OMEGA,BC,FEMSYS,COEFF,HMESH,U,UX) compute 
%   the L2 and the H1 error with respect to the exact solution, 
%   precised in the string U, with derivative in the string UX. 
%   The result is stored in the vector ERROR.
%
%    [UH,COORD,ERROR,FLAG] = FEM1D(OMEGA,BC,FEMSYS,COEFF,HMESH,ADAPT,U,UX) 
%   is the most general form of this function. 
%
%   About the boundary conditions, see IMPOSEBC.
%
 
%   Fausto Saleri, Alessandro Veneziani  2000.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:48 $
%    from 2.0 to 2.2: Robin conditions

% Check for an acceptable number of input arguments
flag = 0;

if nargin < 5
  es=sprintf(['Not enough input arguments.']);
  errore(es); flag = -1;
  return;
end
if nargin >= 9
  es=sprintf(['Too many input arguments.']);
  errore(es); flag = -1;
  return;
end

fem = femsys(1)
types = femsys(2)

switch nargin
case 6,
   adapt = arg1;
case 7,
   u = arg1; ux = arg2; adapt = 0;
case 8,
   adapt = arg1; u = arg2; ux = arg3;
otherwise
   adapt = 0;
end
adapt
     
if ischar(hmesh) == 0   
   [coord,nvert,flag] = uniform_mesh (Omega,hmesh,flag);
else   
   [coord,nvert,flag] = nonunif_mesh (Omega,hmesh,flag);
end

ne = nvert-1;
[dof,ndof,nln,flag] = create_dof (fem,coord,nvert,flag);
[list_el,flag] = create_list (fem,ne,flag);
[matfem, bfem] = assembla(fem,ne,ndof,nln,dof,list_el,coeff,'c ');
[uh]  = imposebc(matfem,bfem,ndof,bc,types);
coord = dof;

if nargin >= 7
   [e]=h1errore(fem, ne, ndof, nln, dof, list_el, uh, u, ux);   
  fprintf(' Error in norm   H1 %g, in norm  L2, %g \n',e(1),e(3));
  out1 = e; out2 = 0;
end

if adapt == 1 & fem == 1
   toll = 1000
   adapt
   fem
   [snorm2]=stimas2(ne,coord,uh);
   newcoord=[];
   for ie = 1:ne
      if snorm2(ie) > (toll/ne)^2
         newcoord=[newcoord,[coord(ie),(coord(ie)+coord(ie+1))/2]];
      else
         newcoord=[newcoord,[coord(ie)]];
      end 
   end
   newcoord=[newcoord,coord(ne+1)]
   %
   % Parte da mettere in una subroutine
   %
   nvert=length(newcoord)
   coord=newcoord;
   ne = nvert-1; flag=0;
   [dof,ndof,nln,flag] = create_dof (fem,coord,nvert,flag);
   [list_el,flag] = create_list (fem,ne,flag);
   [matfem, bfem]=assembla(fem,ne,ndof,nln,dof,list_el,coeff,f);
 if length(bc)==2
   %
   % Dirichlet boundary conditions
   bfem = bfem - matfem (:,1)*bc(1) - matfem (:,ndof)*bc(2);
   bfem = bfem([2:ndof-1]);
   matfem = matfem([2:ndof-1],[2:ndof-1]);
   ufem = matfem\bfem;
   uh = [bc(1); ufem; bc(2)]
 elseif bc(3) < 0
   %
   % Dirichlet type in the last node  <=> Robin/Neumann in the first
   bfem = bfem - matfem (:,ndof)*bc(2);
   bfem = bfem([1:ndof-1]);
   bfem(1) = bfem(1) - bc(1);
   matfem = matfem([1:ndof-1],[1:ndof-1]);
   matfem(1,1) = matfem(1,1) - bc(4);
   ufem = matfem\bfem;
   uh = [ufem;bc(2)];
 elseif bc(3) > 0
   %
   % Dirichlet type in the first node <=> Robin/Neumann in the last
   bfem = bfem - matfem (:,1)*bc(1);
   bfem = bfem(ndof) + bc(2);
   bfem = bfem([2:ndof]);
   matfem = matfem([2:ndof],[2:ndof]);
   matfem(ndof,ndof) = matfem(ndof,ndof) + bc(5);
   ufem = matfem\bfem;
   uh = [bc(1); ufem];
 else
   %
   % Robin (hopefully!) on both the endpoints
   bfem(1) = bfem(1) - bc(1);
   bfem(ndof)= bfem(ndof) + bc(2);
   matfem(1,1) = matfem(1,1) - bc(4);
   matfem(ndof,ndof) = matfem(ndof,ndof) + bc(5);
   uh = matfem\bfem;
 end
 coord = dof;
 if nargin == 9    
   [e]=h1errore(fem, ne, ndof, nln, dof, list_el, uh, u, ux);   
   fprintf(' Error  in norm  H1 %g, in norm  L2, %g \n',e(1),e(3)); 
     out1 = e; out2 = 0;
 end
%
% Fine parte da subroutine
%
end
%
% Ci son da prevedere anche le basi gerarchiche
%

if flag < 0
   uh=[]; coord=[]; dof=[]; ah=[];
end

return

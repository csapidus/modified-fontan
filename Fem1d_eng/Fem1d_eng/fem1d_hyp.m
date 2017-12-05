function [soluzione,coord,out1,out2,out3]=fem1d_hyp(Time,Omega,bc,value_bc,coef_rob,u_0,femsys,coeff,hmesh,forma,theta,DeltaT,...
                                             flag_const,infl,coeff_supp,arg1,arg2,arg3)
%   
%FEM1D        One-dimensional finite element solver for linear 
%             parabolic problem with non constant coefficients.
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
%    FEMSYS(3) is the local Peclet for the transport phenomena.
%    FEMSYS(3) is the local Peclet for the reactive phenomena.
%
%    HMESH is a macro for the grid stepsize. For uniform grid
%   HMESH is a positive number. 
%    FORMA is an integer equal to 1 for conservative form of the
%   equation, 0 otherwise (see also ELLITTICO).
%   In the output variable COORD is stored the computed grid.
%
%    [UH,COORD] = FEM1D(OMEGA,BC,FEMSYS,COEFF,HMESH,FORMA,ADAPT) realize
%   a simple grid adaption in order to reduce the H1-error tunder
%   the positive tolerance ADAPT.
%
%    [UH,COORD,ERROR] = FEM1D(OMEGA,BC,FEMSYS,COEFF,HMESH,FORMA,U,UX) compute 
%   the L2 and the H1 error with respect to the exact solution, 
%   precised in the string U, with derivative in the string UX. 
%   The result is stored in the vector ERROR.
%
%    [UH,COORD,ERROR,FLAG] = FEM1D(OMEGA,BC,FEMSYS,COEFF,HMESH,FORMA,ADAPT,U,UX) 
%   is the most general form of this function. 
%
%   About the boundary conditions, see IMPOSEBC.
%
 
%   Fausto Saleri, Alessandro Veneziani  2000.

% Check for an acceptable number of input arguments
flag = 0;
nargin

if nargin < 15
  disp('Not enough input arguments.');
  %errore(es); 
  flag = -1;
  return;
end
if nargin >= 19
  disp('Too many input arguments.');
  %errore(es); 
  flag = -1;
  return;
end

fem = femsys(1);
types = femsys(2);
%
% tam = time advancing method => il default e' Backward Euler (numero 2)
%
if length(femsys)>2
   tam=femsys(3);
else 
   tam=2;
end;   
if length(femsys)>3
   lump=femsys(4);
else 
   lump=0;
end;   

switch nargin
case 16,
   adapt = arg1;
case 17,
   u = arg1; ux = arg2; adapt = 0;
case 18,
   adapt = arg1; u = arg2; ux = arg3;
otherwise
   adapt = 0;
end

if ischar(hmesh) == 0   
   [coord,nvert,flag] = uniform_mesh (Omega,hmesh,flag);
else   
   [coord,nvert,flag] = nonunif_mesh (Omega,hmesh,flag);
end

ne = nvert-1;
[dof,ndof,nln,flag] = create_dof (fem,coord,nvert,flag);
[list_el,flag] = create_list (fem,ne,flag);
%
%Matrice di massa (con eventuale lumping)
%
if lump==0, massa = crea_massa(fem,ne,ndof,nln,dof,list_el);
else 
   massa=speye(ndof,ndof);
   for i=2:ndof-1, massa(i,i)=1/2*(dof(i+1)-dof(i-1)); end;
   massa(1,1)=dof(2)-dof(1);
   massa(ndof,ndof)=dof(ndof)-dof(ndof-1);
end;
% inizializzazioni
%
t=Time(1);
i_curr=1;
x=dof;
soluzione(i_curr,:)= eval(u_0);
u_old=soluzione(1,:)';
if forma == 1  % forma conservativa
   [matfem, bfem] = assembla_hyp(fem,ne,ndof,nln,dof,list_el,coeff,coeff_supp,'c ',tam,DeltaT,t);
else % forma non conservativa
   [matfem, bfem] = assembla_hyp(fem,ne,ndof,nln,dof,list_el,coeff,coeff_supp,'nc',tam,DeltaT,t);
end;
%
%%%%%%%% GESTIONE DELLE BC
%%
[nu,beta,sigma,f,flag]=dataextr(coeff); % per come sono strutturate le cose, nu=0
if infl==0 
  if forma==1 
     x=dof(ndof);
     beval=eval(beta);
     matfem(ndof,ndof)=matfem(ndof,ndof)+beval;
   end;
   bc_eff=eval(value_bc);
   bfem = bfem - matfem (:,1)*bc_eff(1);
   bold = bfem([2:ndof]);
   matold = matfem([2:ndof],[2:ndof]);
elseif infl==1
   if forma==1 
     x=dof(1);
     beval=eval(beta);
     matfem(1,1)=matfem(1,1)-beval;
   end;
   bc_eff=eval(value_bc);
   bfem = bfem - matfem (:,ndof)*bc_eff;
   bold = bfem([1:ndof-1]);   
   matold = matfem([1:ndof-1],[1:ndof-1]);
elseif infl==2
   bc_eff(1)=eval(value_bc(1,:));
  bc_eff(2)=eval(value_bc(2,:));
  bfem = bfem - matfem (:,1)*bc_eff(1) - matfem(:,ndof)*bc_eff(2);
  bold = bfem([2:ndof-1]);
  matold = matfem([2:ndof-1],[2:ndof-1]);
end;  
%if bc(1)==1 
% if bc(2)>1  %Dirichlet a sx, Neumann/Robin a dx
%  bc_eff(4)=eval(coef_rob(1,:));
%  bc_eff(5)=eval(coef_rob(2,:));
%  bc_eff(3)=1;
%  bfem = bfem(ndof) + bc_eff(2);
%  bfem = bfem - matfem (:,1)*bc_eff(1);
%  bold = bfem([2:ndof]);
%  matfem(ndof,ndof) = matfem(ndof,ndof) + bc_eff(5);
%  matold = matfem([2:ndof],[2:ndof]);
% else % Dirichlet a sx e dx
%  bfem = bfem - matfem (:,1)*bc_eff(1) - matfem(:,ndof)*bc_eff(2);
%  bold = bfem([2:ndof-1]);
%  matold = matfem([2:ndof-1],[2:ndof-1]);
% end;
%else
% if bc(2)>1  %Robin a sx e dx
%  bc_eff(4)=eval(value_bc(1,:));
%  bc_eff(5)=eval(value_bc(2,:));
%  bc_eff(3)=-1;
%  bfem(1) = bfem(1) + bc_eff(1);
%  bfem(ndof) = bfem(ndof) + bc_eff(2);
%  matfem(1,1) = matfem(1,1) + bc_eff(4);
%  matfem(ndof,ndof) = matfem(ndof,ndof) + bc_eff(5);
%  bold=bfem; matold=matfem;
% else %Dirichlet a dx, Neumann/Robin a sx
%  bc_eff(4)=eval(value_bc(1,:));
%  bc_eff(5)=eval(value_bc(2,:));
%  bc_eff(3)=0
%  bfem = bfem - matfem (:,ndof)*bc_eff(2);
%  bfem(1) = bfem(1) + bc_eff(1);
%  bold = bfem([1:ndof-1]);
%  matfem(1,1) = matfem(1,1) + bc_eff(4);
%  matold = matfem([1:ndof-1],[1:ndof-1]);
% end;
%end;   
%
% ciclo temporale 
%
for t=Time(1)+DeltaT:DeltaT:Time(2)
%
disp('Time:'),t
% operazioni da fare nel caso NON costante (o comunque la prima volta)
% assemblaggio matrici pb spaziale   
if flag_const<=0 | i_curr==1
if forma == 1  % forma conservativa
   [matfem, bfem] = assembla_hyp(fem,ne,ndof,nln,dof,list_el,coeff,coeff_supp,'c ',tam,DeltaT,t);
else % forma non conservativa
   [matfem, bfem] = assembla_hyp(fem,ne,ndof,nln,dof,list_el,coeff,coeff_supp,'nc',tam,DeltaT,t);
end;
if infl==0 
  bc_eff=eval(value_bc);
  if forma==1 
     x=dof(ndof);
     beval=eval(beta);
     matfem(ndof,ndof)=matfem(ndof,ndof)+beval;
   end;
elseif infl==1
   bc_eff=eval(value_bc);
   if forma==1 
     x=dof(1);
     beval=eval(beta);
     matfem(1,1)=matfem(1,1)-beval;
   end;
elseif infl==2   
   bc_eff(1)=eval(value_bc(1,:));
   bc_eff(2)=eval(value_bc(2,:));
 end;  
%bc_eff(1)=eval(value_bc(1,:));
%bc_eff(2)=eval(value_bc(2,:));
%if bc(1)==1 
% if bc(2)>1  
%  bc_eff(4)=eval(coef_rob(1,:));
%  bc_eff(5)=eval(coef_rob(2,:));
%  bc_eff(3)=1
% end;
%else
% if bc(2)>1  
%  bc_eff(4)=eval(value_bc(1,:));
%  bc_eff(5)=eval(value_bc(2,:));
%  bc_eff(3)=-1
% else
%  bc_eff(4)=eval(value_bc(1,:));
%  bc_eff(5)=eval(value_bc(2,:));
%  bc_eff(3)=0
% end;
%end;   
%
end; % da qui in avanti è la parte anche a coeff costanti
%attenzione: matred e bred si riferiscono 
% a matrice e termine noto della parte convettiva
%
[uh,matred,bred,stima_cond]  = imposebc_hyp(matfem,bfem,ndof,bc_eff,massa,matold,bold,u_old,...
                                  theta,DeltaT,types);
%
%
coord = dof;
% Analisi dell'errore
%

if nargin >= 16
   x=dof;
   u_eff = eval(u);
   ux_eff = eval(ux);
   [e]=h1errore(fem, ne, ndof, nln, dof, list_el, uh, u, ux,t);   
   errori(i_curr,:)=e;
%  fprintf(' Error in norm H1 %g, in norm L2, %g \n',e(1),e(3));
end 
%
% assegnamenti prima del nuovo time step
%
i_curr=i_curr+1;
soluzione(i_curr,:)=uh';
matold=matred;
bold=bred;
u_old = uh;
clear matred, bred;
end; % ciclo temporale

if adapt == 1 & fem == 1
   toll = 1000
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
   bfem(1) = bfem(1) + bc(1);
   matfem = matfem([1:ndof-1],[1:ndof-1]);
   matfem(1,1) = matfem(1,1) + bc(4);
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
   bfem(1) = bfem(1) + bc(1);
   bfem(ndof)= bfem(ndof) + bc(2);
   matfem(1,1) = matfem(1,1) + bc(4);
   matfem(ndof,ndof) = matfem(ndof,ndof) + bc(5);
   uh = matfem\bfem;
 end
 coord = dof;
 if nargin == 16   
   [e]=h1errore(fem, ne, ndof, nln, dof, list_el, uh, u, ux);   
   fprintf(' Error in norm H1 %g, in norm L2, %g \n',e(1),e(3)); 
     out1 = e; out2 = 0;
 end
%
% Fine parte da subroutine
%
end
if flag < 0
   uh=[]; coord=[]; dof=[]; ah=[];
end

if nargout==4
   out1=stima_cond; out2=matfem;
elseif nargout==5
   out1=errori; out2=stima_cond; out3=matfem;
end;   
return

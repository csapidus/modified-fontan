function [a,b,mu_sx,mu_dx]=assembla(fem,ne,ndof,nln,dof,list_el,coeff,form,stable,t)
%
%ASSEMBLA       Evaluation 
%
%    [A,B]=ASSEMBLA(FEM,NE,NDOF,NLN,DOF,LIST_EL,COEFF) 
%    return the finite element matrix A and the load vector B.
%
 
%   Saleri Fausto, 1999.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:47 $

[nu,beta,sigma,f,flag]=dataextr(coeff);
[nqn, qnodes, w] = xwgll(7)

[phi,dphi] = bases1d (fem);

a = sparse (ndof,ndof); b = zeros(ndof,1);

omega0 = dof(list_el(1,1));
omega1 = dof(list_el(2,ne));
x=omega0;
nu_sx=eval(nu);
beta_sx=eval(beta);
mu1=nu_sx;
x=omega1;
nu_dx=eval(nu);
beta_dx=eval(beta);
mu2=nu_dx;
%%
% per la stabilizzazione visc art
%%
if stable==2
   x=[omega0:1/ne:omega1];
   beta_ev=eval(beta);
   beta_max=norm(beta_ev,inf);
end;   
   
for ie = 1:ne   
   jac= dof (list_el(2,ie)) - dof (list_el(1,ie));  
   [phiq, dphiq] = valuta (fem, nqn, qnodes, phi, dphi);
   if nargin==10  
    [nq, bq, sq, fq] = valutaf (nqn, qnodes, jac, dof(list_el(1,ie)),...
      nu, beta, sigma, f,t);
   else
    [nq, bq, sq, fq] = valutaf (nqn, qnodes, jac, dof(list_el(1,ie)),...
      nu, beta, sigma, f);
   end;      
   if stable ==2 %visc art
      nq=viscart(nq,beta_max,jac);
      mu1=viscart(nu_sx,beta_max,jac);
      mu2=viscart(nu_dx,beta_max,jac);   
   elseif stable==3 %upwind
      nq=upw(nq,bq,jac);
      mu1=upw(nu_sx,beta_sx,jac);
      mu2=upw(nu_dx,beta_dx,jac);
   elseif stable==4 %sg
      nq=sg(nq,bq,jac);
      mu1=sg(nu_sx,beta_sx,jac);
      mu2=sg(nu_dx,beta_dx,jac);
   elseif stable==5 % mass lumping sul termine reattivo per P1 (effettuato fuori)
      sq=0;
   end;
   %% scalatura bd di Neumann/Robin
   if ie==1, mu_sx = mu1;
   elseif ie==ne, mu_dx = mu2;
   end;
   %%
   if form=='nc'
    [aloc,bloc] = locassembla_nc (nln, nqn, nq, sq, bq, ...
                 jac, w, dphiq, phiq, fq);
   else
    [aloc,bloc] = locassembla_c (nln, nqn, nq, sq, bq, ...
                 jac, w, dphiq, phiq, fq);
   end;           
   for iloc = 1:nln
      i = list_el (iloc,ie);
      for jloc = 1:nln
         j = list_el (jloc,ie);
         a (i,j) = a (i,j) + aloc(iloc,jloc);
      end
      b (i) = b (i) + bloc (iloc);
   end
end
if (form=='c ') %handling the conservative form
  beta_sx,beta_dx
  a(1,1)=a(1,1)+beta_sx;
  a(ndof,ndof)=a(ndof,ndof)+beta_dx;
end;  
return
   
      
   
   
   
   
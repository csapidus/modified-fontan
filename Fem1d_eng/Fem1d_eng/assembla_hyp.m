function [a,b]=assembla_hyp(fem,ne,ndof,nln,dof,list_el,coeff,coeff_supp,form,tam,DeltaT,t)
%
%ASSEMBLA_  Evaluation 
%
%    [A,B]=ASSEMBLA(FEM,NE,NDOF,NLN,DOF,LIST_EL,COEFF) 
%    return the finite element matrix A and the load vector B.
%    for a hypervolic problem
 
%   Saleri Fausto, 1999 + Veneziani Alessandro 2000
%   May 2000

[nu,beta,sigma,f,flag]=dataextr(coeff); % per come sono strutturate le cose, nu=0
if tam==4
[nu2,beta2,sigma2,f2,flag2]=dataextr(coeff_supp);
end;
omega0 = dof(list_el(1,1));
omega1 = dof(list_el(2,ne));
x=omega0;
beta_sx=eval(beta);
x=omega1;
beta_dx=eval(beta);

[nqn, qnodes, w] = xwgll(7);

[phi,dphi] = bases1d (fem);

a = sparse (ndof,ndof); b = zeros(ndof,1);

for ie = 1:ne   
   jac= dof (list_el(2,ie)) - dof (list_el(1,ie));  
   [phiq, dphiq] = valuta (fem, nqn, qnodes, phi, dphi);
   if nargin==12  
    [nq, bq, sq, fq] = valutaf (nqn, qnodes, jac, dof(list_el(1,ie)),...
      nu,beta, sigma, f,t);
   else
    [nq, bq, sq, fq] = valutaf (nqn, qnodes, jac, dof(list_el(1,ie)),...
      nu,beta, sigma, f);
   end;      
%   
%
 if tam==3, nq=jac/2*abs(bq); %<<< UPWIND
 elseif tam==4  %<<< LAX-WENDROFF
    [nq2, bq2, sq2, fq2] = valutaf (nqn, qnodes, jac, dof(list_el(1,ie)),...
      nu2,beta2, sigma2, f2,t);
    nq_eff=DeltaT/2*bq.*bq;
    bq_eff=bq+DeltaT/2*bq2;
    sq_eff=sq-DeltaT/2*sq2;
    fq_eff=fq-DeltaT/2*sq.*fq+DeltaT/2*sq2;
    nq=nq_eff;
    bq=bq_eff;
    sq=sq_eff;
    fq=fq_eff;
 end; 
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

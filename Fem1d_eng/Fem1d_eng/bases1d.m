function [phi,dphi,flag] = bases1d (fem)
%BASES1D        Basis function for one-dimensional 
%               finite elements. 
%
%    [PHI,DPHI]=BASES1D(FEM) return in the
%    array PHI and DPHI the basis function and their derivatives.
%    FEM can be equal to -3,-2,-1 for hyerarchical basis,
%    to 3,2,1 for lagrangian basis respectively of degree
%    abs(FEM).
%
 
%   Saleri Fausto, 1999.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:47 $

if nargout == 3
   flag = 0;
end

switch fem
case 1,   
   phi  = ['1-csi';'csi  '];
   dphi = ['-1';' 1'];
case 2,   
   phi = ['2*csi.^2-3*csi+1';'2*csi.^2-csi    ';'4*csi-4*csi.^2  '];
   dphi = ['4*csi-3';'4*csi-1';'4-8*csi'];
case 3,
   phi = ['-0.5*(9*csi.^3-18*csi.^2+11*csi-2)';...
         '0.5*(9*csi.^3-9*csi.^2+2*csi)     ';...
         '4.5*(3*csi.^3-5*csi.^2+2*csi)     ';...
         '-4.5*(3*csi.^3-4*csi.^2+csi)      ';];
   dphi = ['-0.5*(27*csi.^2-36*csi+11)';...
           '0.5*(27*csi.^2-18*csi+2)  ';...
           '4.5*(9*csi.^2-10*csi+2)   ';...
           '-4.5*(9*csi.^2-8*csi+1)   '];
case -1,   
   phi  = ['1-csi';'csi  '];
   dphi = ['-1';' 1'];
case -2
   phi = ['1-csi       ';'csi         ';'(1-csi).*csi'];
   dphi = ['-1     ';'1      ';'1-2*csi'];
case -3
   phi = ['(1-csi)                ';'         csi           ';'(1-csi).*csi           ';'(1-csi).*csi.*(1-2*csi)'];
   dphi = ['-1                         ';' 1                         ';' 1-2*csi                   ';'(1-2*csi).^2-2*csi.*(1-csi)'];
otherwise
   es=sprintf(['Finite element not avaible.']);
   errore(es); flag = -1;
end

return
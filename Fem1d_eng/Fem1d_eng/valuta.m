function [phiq,dphiq,flag] = valuta (fem, nqn, qnodes, phi, dphi)
%
%VALUTA         Evaluation of the basis functions and their
%               first derivatives in the quadrature nodes. 
%
%    [PHIQ,DPHIQ]=VALUTA(FEM,NQN,QNODES,PHI,DPHI) 
%    return in the array PHIQ and DPHIQ the evaluation of
%    the basis function and their derivatives in the quadrature
%    nodes QNODES.
%
 
%   Saleri Fausto, 1999.

if nargout == 3
   flag = 0;
end

phiq = [];  dphiq = [];
switch abs(fem)
   case 1,      
      for k = 1:nqn        
         csi = qnodes (k);
         phiq = [phiq, [eval(phi(1,:)); eval(phi(2,:))]];
         dphiq = [dphiq, [eval(dphi(1,:));eval(dphi(2,:))]];
      end
   case 2,
      for k = 1:nqn        
         csi = qnodes (k);
         phiq = [phiq, [eval(phi(1,:)); eval(phi(2,:));eval(phi(3,:))]];
         dphiq = [dphiq, [eval(dphi(1,:));eval(dphi(2,:));eval(dphi(3,:))]];         
      end
   case 3,
      for k = 1:nqn        
         csi = qnodes (k);
         phiq = [phiq, [eval(phi(1,:)); eval(phi(2,:));eval(phi(3,:));eval(phi(4,:))]];
         dphiq = [dphiq, [eval(dphi(1,:));eval(dphi(2,:));eval(dphi(3,:));eval(dphi(4,:))]];         
      end
   otherwise
      es=sprintf(['Finite elements avaible are P1, P2 and P3']);     
      disp(es); flag=-1; phiq=[]; dphiq=[]; 
      return;
end

return

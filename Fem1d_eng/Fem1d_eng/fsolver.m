function [x]=fsolver(A,b,types)
%
%FSOLVER        Solve a linear system.
%
%    [X]=FSOLVER(A,B) solve the linear system AX=B
%    with a direct method.
%
%    [X]=FSOLVER(A,B,TYPES) solve the linear system AX=B
%    with a direct method if TYPES is equal 0, with
%    Conjugate Gradients Method if TYPES is equal 1,
%    BiConjugate Gradients Method if TYPES is equal 2,
%    BiConjugate Gradients Stabilized Method if TYPES is equal 3,
%    Conjugate Gradients Squared Method if TYPES is equal 4,
%    Generalized Minimum Residual Method if TYPES is equal 5,
%    Quasi-Minimal Residual Method if TYPES is equal 6.
%
 
%   Saleri Fausto, 1999.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:48 $

if nargin == 2
   types = 0;
end

tol = 1.e-10;
maxit = 4*max(size(A));

switch types
case 0,
   [L,U,P] = lu(A);
   y = L\(P*b);
   x = U\y;
case 1,
   [x,flag,relres,iter] = pcg(A,b,tol,maxit);
   switch flag
   case 0,
      fprintf('PCG converged to the desired tolerance %e \n within %d iterations\n',tol,iter);
   case 1,
      fprintf('PCG iterated %d times but did not converge.\n',maxit);
   case 2,
      fprintf('The system is too ill conditioned.\n');
   case 3,
      fprintf('PCG stagnated.\n');
   case 4,
      fprintf('One of the scalar quantities calculated during \n PCG became too small or too large to continue computing.\n')
   end      
case 2,
   [x,flag,relres,iter] = bicg(A,b,tol,maxit);   
   switch flag
   case 0,
      fprintf('BICG converged to the desired tolerance %e \n within %d iterations\n',tol,iter);
   case 1,
      fprintf('BICG iterated %d times but did not converge.\n',maxit);
   case 2,
      disp('The system is too ill conditioned.\n');
   case 3,
      disp('BICG stagnated.\n');
   case 4,
      disp('One of the scalar quantities calculated during \n BICG became too small or too large to continue computing.\n')
   end      
case 3,
   [x,flag,relres,iter] = bicgstab(A,b,tol,maxit);
   switch flag
   case 0,
      fprintf('BICGSTAB converged to the desired tolerance %e \n within %d iterations\n',tol,iter);
   case 1,
      fprintf('BICGSTAB iterated %d times but did not converge.\n',maxit);
   case 2,
      disp('The system is too ill conditioned.\n');
   case 3,
      disp('BICGSTAB stagnated.\n');
   case 4,
      disp('One of the scalar quantities calculated during \n BICGSTAB became too small or too large to continue computing.\n')
   end      
case 4,
   [x,flag,relres,iter] = cgs(A,b,tol,maxit);
   switch flag
   case 0,
      fprintf('CGS converged to the desired tolerance %e \n within %d iterations\n',tol,iter);
   case 1,
      fprintf('CGS iterated %d times but did not converge.\n',maxit);
   case 2,
      disp('The system is too ill conditioned.\n');
   case 3,
      disp('CGS stagnated.\n')
   case 4,
      disp('One of the scalar quantities calculated during \n CGS became too small or too large to continue computing.\n')
   end      
case 5,
   [x,flag,relres,iter] = gmres(A,b,tol,maxit);
   switch flag
   case 0,
      fprintf('GMRES converged to the desired tolerance %e \n within %d iterations\n',tol,iter);
   case 1,
      fprintf('GMRES iterated %d times but did not converge.\n',maxit);
   case 2,
      disp('The system is too ill conditioned.');
   case 3,
      disp('GMRES stagnated.');
   case 4,
      disp('One of the scalar quantities calculated during GMRES became too small or too large to continue computing.')
   end      
case 6,
   [x,flag,relres,iter] = qmr(A,b,tol,maxit);
      switch flag
   case 0,
      fprintf('QMR converged to the desired tolerance %e \n within %d iterations\n',tol,iter);
   case 1,
      fprintf('QMR iterated %d times but did not converge.\n',maxit);
   case 2,
      disp('The system is too ill conditioned.');
   case 3,
      disp('QMR stagnated.');
   case 4,
      disp('One of the scalar quantities calculated during QMR became too small or too large to continue computing.')
   end      
otherwise
   [L,U,P] = lu(A);
   y = L\(P*b);
   x = U\y;
end

return

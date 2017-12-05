function exe_diffox
%
% EXE_DIFFOX
%   execute fem1d program for the slider test case.
%
%   See also FEM1D.

%   Reference: A. Quarteroni, Modellistica Numerica 
%   per Problemi Differenziali, Springer-Italia, Milano
%   (2000).

%   Fausto Saleri, Alessandro Veneziani
%   $Revision: 1.1.1.1 $  $Date: 2001/03/09 08:22:48 $
global CHOICE FORMA OMEGA HMESH FEMSYS BCVALUE
global uh coord
%
% Calcola i parametri per Fem1d
Omega(1) = 0;
Omega(2) = str2num(get(OMEGA(2),'String'));
cext=str2num(get(BCVALUE(1,2),'String'));
bc = [0 cext -1 0];
hmesh = get(HMESH,'String');
femsys(1) = get(FEMSYS,'Value');
femsys(2) = 0;
kappa = get(FORMA(1,3),'String');
mu = get(FORMA(1,1),'String');
mu = ['x^2*',mu];
gamma = ['x^2*',kappa];
coeff = ['mu = ',mu,', beta = 0, gamma =',gamma,', f = 0'];
forma = get(CHOICE(1),'Value');
[uh,coord]=fem1d_ell(Omega,bc,femsys,coeff,hmesh,forma);
return

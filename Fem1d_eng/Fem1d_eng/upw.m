function visc=upw(nq,bq,h);
%
% Viscosit� Upwind
%
global FORMA CHOICE

visc=nq+0.5*h*abs(bq);
return;
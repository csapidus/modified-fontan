function visc=upw(nq,bq,h);
%
% Viscosità Upwind
%
global FORMA CHOICE

visc=nq+0.5*h*abs(bq);
return;
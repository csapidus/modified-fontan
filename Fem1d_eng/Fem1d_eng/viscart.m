function visc=viscart(nq,beta_max,h);
%
% Viscosità Upwind
%
global FORMA CHOICE

visc=nq+h*beta_max;
return;
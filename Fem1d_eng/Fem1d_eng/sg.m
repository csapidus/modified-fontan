function visc=sg(nq,bq,h);
%
% Viscosità Scharfetter-Gummel
%

phi=bern(h*abs(bq)./nq);
visc=0.5*h*abs(bq)+nq.*phi;
return;

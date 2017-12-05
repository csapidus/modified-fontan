function [p,pd,pm1,pdm1,pm2,pdm2] = jacobf(n,alpha,beta,x) 
% 
% 
% 
apb=alpha+beta; 
p=1; 
pd=0; 
if n == 0 
 return 
end 
polyl = p; 
pderl = pd; 
p = (alpha-beta+(apb+2)*x)*.5; 
pd = .5*(apb+2); 
if n == 1 
 return 
end 
for k = 2:n 
  a1 = 2*k*(k+apb)*(2*k+apb-2); 
  a2 = (2*k+apb-1)*(alpha^2-beta^2); 
  b3 = 2*k+apb-2; 
  a3 = b3*(b3+1)*(b3+2); 
  a4 = 2*(k+alpha-1)*(k+beta-1)*(2*k+apb); 
  polyn = ((a2+a3*x)*p-a4*polyl)/a1; 
  pdern = ((a2+a3*x)*pd-a4*pderl+a3*p)/a1; 
  psave = polyl; 
  pdsave = pderl; 
  polyl = p; 
  p = polyn; 
  pderl = pd; 
  pd = pdern; 
end 
pm1 = polyl; 
pdm1 = pderl; 
pm2 = psave; 
pdm2 = pdsave; 
% 
return 
 

function [xjac] = jacg (np,alpha,beta) 
% 
% 
% 
eps=1.e-12; 
kstop=10; 
n=np-1; 
dth = pi*.5/(n+1); 
for j =1:np 
  if j == 1  
    x = cos((2*(j-1)+1)*dth); 
  else 
    x1 = cos((2*(j-1)+1)*dth); 
    x2 = xlast; 
    x = (x1+x2)*.5; 
  end 
  for k =1:kstop 
    [p,pd,pm1,pdm1,pm2,pdm2] = jacobf(np,alpha,beta,x); 
    recsum=0; 
    jm=j-1; 
    for i =1:jm 
     recsum=recsum+1/(x-xjac(np-i+1)); 
    end 
    delx = -p/(pd-recsum*p); 
    x=x+delx; 
    if  abs(delx) < eps 
     break 
    end 
  end 
  xjac(np-j+1) = x; 
  xlast = x; 
end 
for i =1:np 
  xmin =2; 
  for j = i:np 
    if xjac(j) < xmin 
     xmin = xjac(j); 
     jmin = j; 
    end 
  end 
  if jmin ~= i 
   swap = xjac(i); 
   xjac(i) =  xjac(jmin); 
   xjac(jmin) = swap; 
  end 
end 
% 
return 

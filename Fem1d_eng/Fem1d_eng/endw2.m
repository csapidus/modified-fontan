function [w2] = endw2(n,alpha,beta)                    

%

%

%

apb = alpha+beta;                                               

if n == 0                                                  

   w2 = 0;                                                   

   return                                                         

end

ga = gammaf(alpha+1);

gb = gammaf(beta+2);

gab = gammaf(apb+3);                                                             

f1 = ga*gb/gab*(apb+2)*2^(apb+2)*.5;                           

if n == 1                                                 

  w2 = f1;                                                     

  return                                                         

end                                                            

fint1 = ga*gb/gab*2^(apb+2);

ga = gammaf(alpha+2);

gab = gammaf(apb+4);                           

fint2 = ga*gb/gab*2^(apb+3);                                    

f2 = (2*(alpha+2)*fint1-(apb+4)*fint2)*(apb+3)*.25;                                       

if n == 2                                                 

   w2 = f2;                                                    

   return                                                         

end                                                            

for i=3:n                                                      

  di = i-1;                                               

  abn = alpha+beta+di;                                           

  abnn = abn+di;                                                  

  a1 = -(2*(di+alpha)*(di+beta))/(abn*abnn*(abnn+1));      

  a2 = (2*(alpha-beta))/(abnn*(abnn+2));                   

  a3 = (2*(abn+1))/((abnn+2)*(abnn+1));                

  f3 =  -(a2*f2+a1*f1)/a3;                                      

  f1 = f2;                                                      

  f2 = f3;                                                      

end

w2  = f3;

%                                                       

return   
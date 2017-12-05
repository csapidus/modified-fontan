function res=bern(x)
%
% Funzione per il calcolo della Bernoulli function
%
[m,n]=size(x);
dim=max(m,n);
for i=1:dim
   if x(i)==0
      res(i)=1;
   else
      res(i)=abs(x(i))/(exp(abs(x(i)))-1);
   end;    
 end;  
return

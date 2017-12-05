function [f]=heaviside(x)
f=zeros(size(x));
f(find(x>0.5))=1;
return;





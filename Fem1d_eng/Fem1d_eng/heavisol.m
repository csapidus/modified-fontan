function [u]=heavisol(x)
great=find(x>0.5);
u=1/8*x;	
u(great)=-0.5.*x(great).^2+5/8*x(great)-1/8;
%if x<=0.5 
%   u=1/8*x;
%else		
%   u=-0.5*x^2+5/8*x-1/8;
%end	
return


function [u] =heavider(x)
great=find(x>0.5);
u=1/8*ones(size(x));
u(great)=-x(great)+5/8;
return


function [u]=rod_exact(x)
% [u]=rod_exact(x) 
%exact solution of problem 12.1 of Numerical Mathematics Book
%
u0=10;
sigma=2;
mu=200;
A=4*pi;
p=4*pi;
L=100;
m=sqrt(sigma*p/(mu*A));
fact=u0/cosh(m*L);
u=fact*cosh(m*(L-x));
return;

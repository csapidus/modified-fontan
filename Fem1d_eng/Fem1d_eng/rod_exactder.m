function [u]=rod_exactder(x)
% [u]=rod_exactder(x) 
%derivative of exact solution of problem 12.1 of Numerical Mathematics Book
%
u0=10;
sigma=2;
mu=200;
A=4*pi;
p=4*pi;
L=100;
m=sqrt(sigma*p/(mu*A));
fact=-m*u0/cosh(m*L);
u=fact*sinh(m*(L-x));
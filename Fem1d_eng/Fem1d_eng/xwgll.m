function [nqn,qnodes,w] = xwgll(np)

% 

% Scopo:  Calcolo degli n nodi e pesi della formula

%       di quadratura di Gauss Lobatto Legendre in [-1,1].

%          

% Uso:	 [x,w] = xwgll(np)

%

% Parametri di ingresso

%  np    numero di punti

%

% Parametri in uscita

%  x    nodi della formula di quadratura

%  w    pesi della formula di quadratura

%

% Autori: pg

%

nqn=np;
if np<=1

  disp('The minimum sumber of quadrature nodes (GLL) is  2');

  return

end

alpha=0;

beta=0;

n=np-1; 

nm1=n-1;

if nm1 > 0 

  alpha1=alpha+1;

  beta1=beta+1;

  [x2,w2] = zwgjd(nm1,alpha1,beta1);

end

x(1) = -1;

x(2:n)=x2;

x(np) = 1;

% pesi

%for i = 1:np

% p=pnleg(x(i),n);

% w(i)=2/(n*(n+1))*1/(p^2);

%end  

[p,pd,pm1,pdm1,pm2,pdm2]=jacobf(n,alpha,beta,x(1));

ww=endw1(n,alpha,beta);

w(1)=ww/(2*pd);

w(2:n)=w2./(1-x2.^2);

[p,pd,pm1,pdm1,pm2,pdm2]=jacobf(n,alpha,beta,x(np));

ww=endw2(n,alpha,beta);

w(np)=ww/(2*pd);

%
% normalizzazione fra 0 e 1
qnodes = 1/2*(x+1);
w = w'/2;


return

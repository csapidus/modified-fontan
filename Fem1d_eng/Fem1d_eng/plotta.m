figure(1)
x=[0:0.01:1];
f=1 - 3/2*x + 9/8*x.^2;
f=[f 1]; x=[x 1];
patch(x,f,'red'); patch([0 1 1 0],[-0.1 -0.1 0 0],'red'); axis off
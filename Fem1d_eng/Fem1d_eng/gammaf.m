function [g] = gammaf(x)

%

%

%

g=1;

if x == -.5

 g = -2*sqrt(pi);

elseif x == .5 

 g = sqrt(pi);

elseif x == 1 | x == 2

 g = 1;

elseif x == 1.5

 g = sqrt(pi)*.5;

elseif x == 2.5

 g = .75*sqrt(pi);

elseif x == 3

 g = 2;

elseif x == 3.5

 g = 15/8*sqrt(pi);

elseif x == 4

 g = 6;

elseif x == 5

 g = 24;

elseif x == 6

 g = 120;

end

%

return

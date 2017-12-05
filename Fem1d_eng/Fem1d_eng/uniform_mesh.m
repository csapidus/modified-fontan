function [coord,nvert,flag] = uniform_mesh(Omega,hmesh,flag)
%   
% uniform_mesh    Uniform grid with grid-space equal to hmesh. 
%
%    [coord,nvert] = uniform_mesh(Omega,hmesh) returns in the
%    1 x nvert array the coordinates of the nodes of the mesh
%    for the domain [Omega(1),Omega(2)].
%
 
%   Saleri Fausto, 1999.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:48 $

if flag < 0
   coord=[]; nvert=0;
   return; 
end

if Omega(1) >= Omega(2)
  es=sprintf(['Omega(1) should be less than Omega(2)']);
  disp(es); flag=-1; coord=[]; nvert=0;
  return;
end
   
if hmesh <= eps*1.e03
   es=sprintf(['hmesh should be greater than %0.5g'],eps*1.e03);
   disp(es); flag=-1; coord=[]; nvert=0;
   return;
end
   
coord=[Omega(1)]; last=Omega(1);

while last <= Omega(2)-hmesh
   last=last+hmesh;
   coord=[coord,last];
end

nvert = length(coord);
if Omega(2)-coord(nvert) < hmesh/5
   coord(nvert)=Omega(2);
else
   nvert=nvert+1;
   coord=[coord,Omega(2)];
end

return

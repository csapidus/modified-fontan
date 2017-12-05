function [list_el,flag] = create_list(fem,ne,flag)
%   
%CREATE_LIST    Reference matrix for the finite element grid. 
%
%    [LIST_EL]=CREATE_LIST(FEM,NE,FLAG) stores in the 
%    NPL x NE array LIST_EL the indices of the nodes
%    of any grid element.
 
%   Saleri Fausto, 1999.
%   $Revision: 1.1.1.1 $ $Date: 2001/03/09 08:22:47 $

switch nargin
case 3,   
   if flag < 0, list_el = []; return; end
otherwise,
   switch nargout
   case 2,
      flag = 0;
   end
end

switch abs(fem)
case 1, list_el = [[1:ne];[2:ne+1]];
case 2, i=[1:ne]; list_el = [ 2*i-1; 2*i+1; 2*i];
case 3, i=[1:ne]; list_el = [ 3*i-2; 3*i+1; 3*i-1; 3*i];
otherwise,   
   es=sprintf(['Finite elements available are P1, P2 and P3']);
   disp(es); flag=-1; list_el=[]; return;
end

return
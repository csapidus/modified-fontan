function call_xor(CHI,cosa,come_si,come_no)
%
%CALL_XOR: funzione di servizio per interfaccia grafica di fem1d:
% fa comparire e scomparire parti di interfaccia a seconda dei casi
%
%   See also FEM1D.

%   Fausto Saleri, Alessandro Veneziani
%   $Revision: 1.1.1.1 $  $Date: 2001/03/09 08:22:47 $
[m,n]=size(CHI);
if length(get(CHI(1),cosa))==length(come_si)
 if get(CHI(1),cosa)==come_si
   for i=1:n
      set(CHI(i),cosa,come_no);
   end;
 else
   for i=1:n
      set(CHI(i),cosa,come_si);
   end;
 end;  
else   
   for i=1:n
      set(CHI(i),cosa,come_si);
   end;
end;
return;

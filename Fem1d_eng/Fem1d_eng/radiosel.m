function [selected]=radiosel(child)
%
%RADIO_SEL .
%   RADIO_SEL .
%
%   See also .

%   Fausto Saleri, Alessandro Veneziani
%   $Revision: 1.1.1.1 $  $Date: 2001/03/09 08:22:48 $
n = max(size(child));
selected = 1;
for i = 1:n
   if get(child(i),'Value') == 1
      selected = i+1;
   end
end
return
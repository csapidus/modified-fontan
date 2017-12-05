function [out_vec] = comp_strings(stringa1,stringa2);
%
% Scopo: prendere stringhe di lunghezza diversa, confrontarle
% e aggiungere gli spazi bianchi a quelle corte per
% compattare tutto in un vettore
%
%A. Veneziani, Aprile 2000
%
l1 = length(stringa1)
l2 = length(stringa2)
if l1>l2
   out_vec(1,:) = stringa1;
   for i=1:l1-l2
    stringa2=[stringa2,' '];
   end;
   out_vec(2,:) = stringa2;
elseif l1==l2
   out_vec(1,:) = stringa1;
   out_vec(2,:) = stringa2;
else
   out_vec(2,:) = stringa2;
   for i=1:l2-l1
    stringa1=[stringa1,' '];
   end;
   out_vec(1,:) = stringa1;
end;
return;
function [der_conom_sx,der_conom_dx]=der_co_nom(choice)
%
%
%
%
global FORMA OMEGA
if choice==1
   str1=get(FORMA(1,1),'String');
   str2=get(FORMA(1,2),'String');
   strx=get(OMEGA(1),'String');
   x=str2num(strx);
   nu=eval(str1);
   beta=eval(str2);
   str1=num2str(nu);
   str2=num2str(beta);
   if str1=='1'
      str1=' ';
   end   
   str1=cat(2,'-',str1);
   str2=cat(2,'+',str2);
   str1=cat(2,str1,'u_x');
   str2=cat(2,str2,'u  ');
   if get(FORMA(1,2),'String')=='0'
      str2='  ';
   elseif  get(FORMA(1,2),'String')=='1'
      str2='+u ';
   end
   der_conom_sx=cat(2,str1,str2);
   if length(der_conom_sx)>13
      der_conom_sx='-m u_x + b u';
   end   
   str1=get(FORMA(1,1),'String');
   str2=get(FORMA(1,2),'String');
   strx=get(OMEGA(2),'String');
   x=str2num(strx);
   nu=eval(str1);
   beta=eval(str2);
   str1=num2str(nu);
   str2=num2str(beta);
   if str1=='1'
      str1=' ';
   end   
   str1=cat(2,' ',str1);
   str2=cat(2,'-',str2);
   str1=cat(2,str1,'u_x');
   str2=cat(2,str2,'u  ');
   if get(FORMA(1,2),'String')=='0'
      str2='  ';
   elseif  get(FORMA(1,2),'String')=='1'
      str2='-u ';
   end
   der_conom_dx=cat(2,str1,str2);
   if length(der_conom_dx)>13
      der_conom_dx='m u_x - b u';
   end   
else
   str0=get(FORMA(2,1),'String');
   str1=str0;
   strx=get(OMEGA(1),'String');
   x=str2num(strx);
   nu=eval(str1);
   str1=num2str(nu);
   if str1=='1'
      str1=' ';
   end   
   str1=cat(2,'-',str1);
   str1=cat(2,str1,'u_x');
   if length(str1)>13
      str1='-m u_x';
   end   
   der_conom_sx=str1;
   str1=str0;
   strx=get(OMEGA(2),'String');
   x=str2num(strx);
   nu=eval(str1);
   str1=num2str(nu);
   if str1=='1'
      str1=' ';
   end   
   str1=cat(2,' ',str1);
   str1=cat(2,str1,'u_x');
   if length(str1)>13
      str1=' m u_x';
   end   
   der_conom_dx=str1;
end
return
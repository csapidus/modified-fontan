function delta=delta_t_set;
%
% selezione del delta_t in base alla CFL
%
global CHOICE FORMA OMEGA METODO HMESH LUMPING

hmesh=get(HMESH,'String')
tam=get(METODO,'Value')

forma=get(CHOICE(1),'Value');
if forma==0, forma=2; end;
find_x = find(hmesh=='x');
if isempty(find_x) == 1
  h = str2num(hmesh);     % hmesh costante
  Omega(1)=str2num(get(OMEGA(1),'String'));
  Omega(2)=str2num(get(OMEGA(2),'String'));
  x=[Omega(1):h:Omega(2)]; 
else
   Omega(1)=str2num(get(OMEGA(1),'String'));
   Omega(2)=str2num(get(OMEGA(2),'String'));
   [x,nvert,flag] = nonunif_mesh (Omega,hmesh,0);
   [m,n]=size(x);
   n=max(m,n);
   h=min(x(2:n)-x(1:n-1));
end
str_coef=[get(FORMA(forma,1),'String')];
num_coef=eval(str_coef);
coef=max(abs(num_coef));
%coef=str2num(get(FORMA(forma,1),'String'));
if (tam==4) & get(LUMPING,'Value')==0
 delta=h/abs(coef)*1/sqrt(3);
elseif (tam==3) & get(LUMPING,'Value')==0
 delta=h/abs(coef)*1/3;
elseif (tam==1) 
   delta=h;
elseif (tam==2)
   delta=h; 
else
   get(LUMPING,'Value')
   delta=h/abs(coef)
end;   

return;

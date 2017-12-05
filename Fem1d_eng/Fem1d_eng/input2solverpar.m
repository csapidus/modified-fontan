function input2solverpar
%
% Funzione di lancio del problema parabolico.
% Modifica del caso ellittico (radicale dal punto di vista informatico)
%
% A. Veneziani - Aprile 2000
%
global OMEGA BC_SX BC_DX BCVALUE HMESH GERBASE FEMSYS
global CHOICE FORMA SOLVER FATHER INFO METODO MODUS
global USCITA AN_ERR AQUAMARINE
global INITC THETA TIME DELTAT
global coord uh mat_fem coord_old uh_old stringa_uex stringa_uxex
%
% Estremi del dominio di calcolo
Omega(1) = str2num(get(OMEGA(1),'String'));
Omega(2) = str2num(get(OMEGA(2),'String'));
Time(1) = str2num(get(TIME(1),'String'));
Time(2) = str2num(get(TIME(2),'String'));
%
% Condizioni iniziali 
%
u0_str=get(INITC,'String');
u0_str=[u0_str,' + 0.*x'];
%
% Condizioni al contorno 
%
% bc è una tabella: 1=Dirichlet; 2=Neumann; 3=Robin
% value_bc = STRINGHE con le b.c. imposte (in generale funzione del tempo
% coef_rob = STRINGHE con gli eventuali coef di Robin
bc = [0; 0];
coef_rob=['0*t';'0*t'];
if  get(BC_SX(1),'Value')== 1
  bc(1) = 1; 
  value_sx = [get(BCVALUE(1,1),'String'),' + 0*t'];
elseif get(BC_SX(2),'Value')== 1
   bc(1) = 2;
   value_sx = [get(BCVALUE(2,1),'String'),' + 0*t'];
else
  bc(1) = 3;   
  alpha_sx = [get(BCVALUE(3,1),'String'),' + 0*t'];
  value_sx = [get(BCVALUE(4,1),'String'),' + 0*t'];
end;   
if  get(BC_DX(1),'Value')== 1
  bc(2) = 1; 
  value_dx = [get(BCVALUE(1,2),'String'),' + 0*t'];
elseif get(BC_DX(2),'Value')== 1
   bc(2) = 2;
   value_dx = [get(BCVALUE(2,2),'String'),' + 0*t'];
else
  bc(2) = 3;   
  alpha_dx = [get(BCVALUE(3,2),'String'),' + 0*t'];
  value_dx = [get(BCVALUE(4,2),'String'),' + 0*t'];
end;   
value_bc = comp_strings(value_sx,value_dx);
if bc(1)==3 & bc(2)==3
   coef_rob = comp_strings(alpha_sx,alpha_dx);
elseif bc(1)==3
   coef_rob= comp_strings(alpha_sx,'0*t');
elseif bc(2)==3
   coef_rob=comp_strings('0*t',alpha_dx);
end;
%
% fine della definizione del problema:
% ora si specifica il metodo
% mesh
hmesh = get(HMESH,'String');
find_x = find(hmesh=='x');
if isempty(find_x) == 1
  hmesh = str2num(hmesh);       % hmesh costante
end
%
% specifiche del metodo di avanzamento in tempo
%
DeltaT = str2num(get(DELTAT,'string'));
theta = str2num(get(THETA,'string'));
%
%tipo di base
%
if get(GERBASE,'Value') == 0
  femsys(1) = get(FEMSYS,'Value');      % base lagrangiana
else
  femsys(1) = -get(FEMSYS,'Value');     % base gerarchica
end
[Pe_l_t,Pe_l_r]=comp_Peclet(get(MODUS,'Value'));
%
if get(CHOICE(1),'Value') == 1
  coeff = ['mu = ',get(FORMA(1,1),'String'),...
           ', beta = ',get(FORMA(1,2),'String'),...
           ', gamma = ',get(FORMA(1,3),'String'),...
           ', effe = ',get(FORMA(1,4),'String')];
else
  coeff = ['mu = ',get(FORMA(2,1),'String'),...
           ', beta = ',get(FORMA(2,2),'String'),...
           ', gamma = ',get(FORMA(2,3),'String'),...
           ', effe = ',get(FORMA(2,4),'String')];
end
femsys(2) = get(SOLVER(1),'Value')-1;
% stabilizzazione
femsys(3) = get(METODO,'Value');
%
% forma contiene info sulla forma e sul fatto che i coefficienti siano costanti
%
forma(1) = get(CHOICE(1),'Value');
flag_const = get(CHOICE(3),'Value');

%
% post-processing
%
analisi_errore = get(AN_ERR(1),'Value');
condizionam = get(USCITA(1),'Value');
patter = get(USCITA(3),'Value');

%%%%%%%%%%%%%%%%
if analisi_errore==0
%%%%%%% soluzione effettiva   
  [uh,coord,stima_cond,mat_fem]=fem1d_par(Time,Omega,bc,value_bc,coef_rob,u0_str,femsys,coeff,hmesh,forma,theta,DeltaT,...
     flag_const);
 if condizionam==1
    set(USCITA(2),'String',num2str(stima_cond));
 end;   
 [figurePos]=pos_figure;
 figurePos(1)=figurePos(1)+50;
FATHER(5)=figure('NumberTitle','off','Color',AQUAMARINE,...
      'Units','points',...
      'MenuBar','none','Name',...
      'Numerical Results',...
      'Resize','off','Position',figurePos);
 % Info
INFO=uicontrol('Units','normalized',...
   'Position',[0.4 0.05 0.40 0.08],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','Press   `Enter` to continue   ',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','normal','FontSize',[8],'Visible','off');
subplot('Position',[0.075 0.25 0.875 0.42]);
 if patter==1
    spy(mat_fem);
    set(INFO,'Visible','on');
    pause
    set(INFO,'Visible','off');
 end;  
 if femsys(1)==-2 
    [m,n]=size(uh)
    uh(:,2:2:m-1) = 0.5*(uh(:,1:2:m-2)+uh(:,3:2:m))+0.25*uh(:,2:2:m-1);
 elseif femsys(1)==-3
    [m,n]=size(uh)
    uh_plot=uh;
    uh_plot(2:3:m-2) = (2*uh(:,1:3:m-3)+2/9*uh(:,3:3:m-1)+2/3*uh(:,2:3:m-2)+uh(:,4:3:m))/3;
    uh_plot(3:3:m-1) = (uh(:,1:3:m-3)-2/9*uh(:,3:3:m-1)+2/3*uh(:,2:3:m-2)+2*uh(:,4:3:m))/3;
    uh=uh_plot;   
 end;   
 tempo=[Time(1):DeltaT:Time(2)];
 if get(AN_ERR(10),'Value')==1
   i=0; 
   for t=Time(1):DeltaT:Time(2)
    i=i+1;   
    PLOT=plot(coord,uh(i,:));
    titolo=['Time  =',num2str(t)];
    T=title(titolo);
    pause(1);
   end;   
 else   
    PLOT=mesh(coord,tempo,uh);
    TITLE=title('Computed Solution  '); 
    XL=xlabel('x');
    YL=ylabel('t');
 end;    
 uicontrol('Units','normalized',...
   'Position',[0.1 0.1 0.2 0.05],'HorizontalAlignment','left',...
   'Style','checkbox',...
   'String','OK','FontAngle','oblique',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[14],'CallBack',...
   'coord_old=coord;uh_old=uh;for i=length(FATHER):-1:3,close(FATHER(i)); end;');
else
 stringa_uex=get(AN_ERR(3),'String');
 stringa_uxex=get(AN_ERR(5),'String');
%%%% soluzione effettiva
 [uh,coord,errori,stima_cond,mat_fem]=fem1d_par(Time,Omega,bc,value_bc,coef_rob,u0_str,femsys,coeff,hmesh,forma,theta,DeltaT,...
     flag_const,stringa_uex,stringa_uxex);
 if condizionam==1
    set(USCITA(2),'String',num2str(stima_cond));
 end;   
%%%%% 
e(3) = norm(errori(:,3),inf);
% calcolo della norma l^2(H^1)
e(1) = 0;num_ts=max(size(errori(:,1)));
for i=1:num_ts, e(1)=e(1)+DeltaT*errori(i,1)^2; end;
%e(1)=e(1)+DeltaT*(errori(1,1)/2+errori(num_ts,1)/2);
e(1)=sqrt(e(1));
%
set(AN_ERR(7),'String',num2str(e(3)));
 set(AN_ERR(9),'String',num2str(e(1)));
 if condizionam==1  
  c2=cond(full(mat_fem),2);
  set(USCITA(2),'String',num2str(c2));
 end;   
 [figurePos]=pos_figure;
 figurePos(1)=figurePos(1)+50;
 FATHER(5)=figure('NumberTitle','off','Color',AQUAMARINE,...
      'Units','points',...
      'MenuBar','none','Name',...
      'Numerical Result ',...
      'Resize','off','Position',figurePos);
   % Info
INFO=uicontrol('Units','normalized',...
   'Position',[0.4 0.05 0.40 0.08],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','Press  `Enter` to continue   ',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','normal','FontSize',[8],'Visible','off');
 subplot('Position',[0.075 0.25 0.875 0.42]);
 if patter==1
   spy(mat_fem);
   set(INFO,'Visible','on');
   pause
   set(INFO,'Visible','off');
end;  
  %x=[Omega(1):hmesh/1000:Omega(2)];
  %u_esatta=eval(stringa_uex);
 if femsys(1)==-2
    [m,n]=size(uh);
    uh(2:2:m-1) = 0.5*(uh(1:2:m-2)+uh(3:2:m))+0.25*uh(2:2:m-1);
elseif femsys(1)==-3
    [m,n]=size(uh)
    coord,uh
    uh_plot=uh;
    uh_plot(2:3:m-2) = (2*uh(1:3:m-3)+2/9*uh(3:3:m-1)+2/3*uh(2:3:m-2)+uh(4:3:m))/3;
    uh_plot(3:3:m-1) = (uh(1:3:m-3)-2/9*uh(3:3:m-1)+2/3*uh(2:3:m-2)+2*uh(4:3:m))/3;
    uh=uh_plot;   
 end;   
  tempo=[Time(1):DeltaT:Time(2)];
  x=[Omega(1):hmesh/1000:Omega(2)];
  if get(AN_ERR(10),'Value')==1
   i=0; 
   for t=Time(1):DeltaT:Time(2)
    stringa_uex;  
    u_esatta=eval(stringa_uex);
    i=i+1;   
    PLOT=plot(coord,uh(i,:),x,u_esatta,'g');
    L=legend('Computed Solution  ','Exact Solution  ');     
    titolo=['Time  ='  num2str(t)];
    T=title(titolo);
    pause(0.1);
   end;   
 else   
    PLOT=mesh(coord,tempo,uh);
    TITLE=title('Exact Solution     '); 
    XL=xlabel('x');
    YL=ylabel('t');
 end;    
 uicontrol('Units','normalized',...
   'Position',[0.1 0.1 0.2 0.05],'HorizontalAlignment','left',...
   'Style','checkbox',...
   'String','OK','FontAngle','oblique',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[14],'CallBack',...
   'coord_old=coord;uh_old=uh;for i=length(FATHER):-1:3,close(FATHER(i)); end;');
end 
return

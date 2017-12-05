function input2solver


global OMEGA BC_SX BC_DX BCVALUE HMESH GERBASE FEMSYS
global CHOICE FORMA SOLVER FATHER INFO METODO MODUS
global USCITA AN_ERR AQUAMARINE
global coord uh mat_fem coord_old uh_old stringa_uex stringa_uxex
%
% Estremi del dominio di calcolo
Omega(1) = str2num(get(OMEGA(1),'String'));
Omega(2) = str2num(get(OMEGA(2),'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Condizioni al contorno 
bc = [0; 0];
coef_rob = [0; 0];
sinistra = get(BC_SX(1),'Value');
if sinistra == 1
  bc(1) = str2num(get(BCVALUE(1,1),'String'));    % Dirichlet a sinistra
else
  sinistra = get(BC_SX(2),'Value');
  bc(3) = -1;
  if  sinistra == 1            
    bc(1) = str2num(get(BCVALUE(2,1),'String'));  % Neumann  a sinistra
    bc(4) = 0;
  else
     bc(1) = str2num(get(BCVALUE(4,1),'String'));  % Robin  a sinistra
     bc(4) = str2num(get(BCVALUE(3,1),'String')); 
  end
end
destra = get(BC_DX(1),'Value');
if destra == 1
  bc(2) = str2num(get(BCVALUE(1,2),'String'));    % Dirichlet a destra
else
  destra = get(BC_DX(2),'Value');
  if length(bc)>2, bc(3) = 0; else, bc(3) = 1; end
   if  destra == 1
     bc(2) = str2num(get(BCVALUE(2,2),'String'));  % Neumann a destra
     bc(5) = 0;    
   else
     bc(5) = str2num(get(BCVALUE(3,2),'String')); 
     bc(2) = str2num(get(BCVALUE(4,2),'String'));  % Robin a destra
  end
end
hmesh = get(HMESH,'String');
find_x = find(hmesh=='x');
if isempty(find_x) == 1
  hmesh = str2num(hmesh);       % hmesh costante
end
if get(GERBASE,'Value') == 0
  femsys(1) = get(FEMSYS,'Value');      % base lagrangiana
else
  femsys(1) = -get(FEMSYS,'Value');     % base gerarchica
end
[Pe_l_t,Pe_l_r]=comp_Peclet(get(MODUS,'Value'));
%
if get(CHOICE(1),'Value') == 1
  coeff = ['mu =     ',get(FORMA(1,1),'String'),...
           ', beta = ',get(FORMA(1,2),'String'),...
           ', gamma= ',get(FORMA(1,3),'String'),...
           ', effe = ',get(FORMA(1,4),'String')];
else
  coeff = ['mu =     ',get(FORMA(2,1),'String'),...
           ', beta = ',get(FORMA(2,2),'String'),...
           ', gamma= ',get(FORMA(2,3),'String'),...
           ', effe = ',get(FORMA(2,4),'String')];
end
femsys(2) = get(SOLVER(1),'Value')-1;
% stabilizzazione
femsys(3) = get(METODO,'Value');
forma = get(CHOICE(1),'Value');
analisi_errore = get(AN_ERR(1),'Value');
condizionam = get(USCITA(1),'Value');
patter = get(USCITA(3),'Value');
if analisi_errore==0
 if (condizionam==1)|(patter==1)
   [uh,coord,mat_fem]=fem1d_ell(Omega,bc,femsys,coeff,hmesh,forma);
 else
   [uh,coord]=fem1d_ell(Omega,bc,femsys,coeff,hmesh,forma);
 end
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
   'String','Type ENTER to continue',...
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
    uh(2:2:m-1) = 0.5*(uh(1:2:m-2)+uh(3:2:m))+0.25*uh(2:2:m-1);
 elseif femsys(1)==-3
    [m,n]=size(uh)
    coord,uh
    uh_plot=uh;
    uh_plot(2:3:m-2) = (2*uh(1:3:m-3)+2/9*uh(3:3:m-1)+2/3*uh(2:3:m-2)+uh(4:3:m))/3;
    uh_plot(3:3:m-1) = (uh(1:3:m-3)-2/9*uh(3:3:m-1)+2/3*uh(2:3:m-2)+2*uh(4:3:m))/3;
    uh=uh_plot;   
 end;   
 if get(AN_ERR(10),'Value')==1
    PLOT=plot(coord,uh,coord_old,uh_old,'m'); 
    L=legend('Numerical Solution ','Computed Solution   ');    
 else   
    PLOT=plot(coord,uh);
    TITLE=title('Computed Solution  '); 
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
 if (condizionam==1)|(patter==1)
  [uh,coord,e,mat_fem]=fem1d_ell(Omega,bc,femsys,coeff,hmesh,forma,stringa_uex,stringa_uxex);
 else
  [uh,coord,e]=fem1d_ell(Omega,bc,femsys,coeff,hmesh,forma,stringa_uex,stringa_uxex);
end
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
      'Risutati Numerici',...
      'Resize','off','Position',figurePos);
   % Info
INFO=uicontrol('Units','normalized',...
   'Position',[0.4 0.05 0.40 0.08],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','premere `Invio` per continuare',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','normal','FontSize',[8],'Visible','off');
 subplot('Position',[0.075 0.25 0.875 0.42]);
 if patter==1
   spy(mat_fem);
   set(INFO,'Visible','on');
   pause
   set(INFO,'Visible','off');
end;  
  x=[Omega(1):hmesh/1000:Omega(2)];
  u_esatta=eval(stringa_uex);
 if femsys(1)==-2 % se con base gerarchica urge una sistemazione in visualizzazione
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
 if get(AN_ERR(10),'Value')==1
    PLOT=plot(coord,uh,coord_old,uh_old,'m',x,u_esatta,'g'); 
    L=legend('Computed Solution  ','Previous Solution   ','Exact Solution  ');    
 else   
    PLOT=plot(coord,uh,x,u_esatta,'g');
    L=legend('Computed Solution  ','Exact Solution  '); 
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

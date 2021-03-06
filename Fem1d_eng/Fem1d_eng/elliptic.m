function elliptic
%
%ELLITTICO 1-d finite element code for elliptic problems.
%   ELLITTICO is the graphical interface for the finite
%   element one-dimensional elliptic solver.
%
%   See also FEM1D.

%   Reference: A. Quarteroni, Modellistica Numerica 
%   per Problemi Differenziali, Springer-Italia, Milano
%   (2000).

%   Fausto Saleri, Alessandro Veneziani
%   $Revision: 1.1.1.1 $  $Date: 2001/03/09 08:22:48 $
global AQUAMARINE RELDATA CHOICE BC_SX BC_DX
global MODUS OMEGA FORMA BCVALUE DER_CON
global FATHER
close all;
%
% Posiziona la figura
[figurePos]=pos_figure;
%
% Crea la figura
FATHER(2)=figure('NumberTitle','off','Color',AQUAMARINE,...
   'Units','points',...
   'MenuBar','none','Name','Fem 1D: elliptic problem',...
   'Resize','off','Position',figurePos);
%
% Costruisce titolo, autori e release
cornici('no ');
uicontrol('Units','normalized',...
   'Position',[0 0.95 1 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',...
   'Formulation and input problem parameters: (constants           )',...
   'BackgroundColor',AQUAMARINE,'FontAngle','oblique',...
   'FontWeight','demi','FontSize',[12]); 
%
% Definizione delle due forme (conservativa e non) del problema.
space = blanks(18);
forma_1 = ['- (',space,'u_x)_x + (',space,'u)_x +',space,'u =',...
      space];
forma_2 = ['- (',space,'u_x)_x  + ',space,'  u_x +',space,'u =',...
      space];
FORMA(1,5)=uicontrol('Units','normalized',...
   'Position',[0.1 0.68 1 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',forma_1,'ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'Visible','off');
FORMA(2,5)=uicontrol('Units','normalized',...
   'Position',[0.1 0.68 1 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',forma_2,'ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'Visible','off');
%
% Bottone di scelta nel caso a coefficienti costanti
CHOICE(3)=uicontrol('Units','normalized',...
   'Position',[0.85 0.96 0.02 0.03],'HorizontalAlignment','left',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[12],...
   'Style','radiobutton','String',' ','Value',1);
%
% Coefficienti della forma conservativa.
FORMA(1,1)=uicontrol('Units','normalized',...
   'Position',[0.13 0.68 0.11 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'Visible','off','callback','if get(FORMA(1,1),stringa)==zero, gest_errori(3), else smista_coef_bc(1); end;');  % viscosita'
FORMA(1,2)=uicontrol('Units','normalized',...
   'Position',[0.38 0.68 0.11 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'Visible','off','callback','smista_coef_bc(1);');  % trasporto
FORMA(1,3)=uicontrol('Units','normalized',...
   'Position',[0.572 0.68 0.11 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'Visible','off');  % reazione
FORMA(1,4)=uicontrol('Units','normalized',...
   'Position',[0.74 0.68 0.11 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'Visible','off');  % termine forzante
%
% Coefficienti della forma non conservativa.
FORMA(2,1)=uicontrol('Units','normalized',...
   'Position',[0.13 0.68 0.11 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1',...
   'BackgroundColor',AQUAMARINE,'ForegroundColor','blue',...
   'FontWeight','demi','FontSize',[10],'Visible','off','callback','if get(FORMA(2,1),stringa)==zero, gest_errori(3), else smista_coef_bc(0); end;');  % viscosita'
FORMA(2,2)=uicontrol('Units','normalized',...
   'Position',[0.38 0.68 0.11 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'Visible','off');  % trasporto
FORMA(2,3)=uicontrol('Units','normalized',...
   'Position',[0.572 0.68 0.11 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'Visible','off');  % reazione
FORMA(2,4)=uicontrol('Units','normalized',...
   'Position',[0.74 0.68 0.11 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'Visible','off');  % termine forzante
%
%
% Bottoni di scelta tra le due forme
CHOICE(1)=uicontrol('Units','normalized',...
   'Position',[0.02 0.86 0.6 0.06],'HorizontalAlignment','left',...
   'BackgroundColor',AQUAMARINE,'Value',0,...
   'FontWeight','demi','FontSize',[12],...
   'Style','radiobutton','String','Conservative Form ','CallBack',...
   'set(CHOICE(1),value,1);set(CHOICE(2),value,0);call_xor(FORMA(1,:),visible,ON,ON);call_xor(FORMA(2,:),visible,OFF,OFF);[dc_sx,dc_dx]=der_co_nom(1);'); 
CHOICE(2)=uicontrol('Units','normalized',...
   'Position',[0.02 0.76 0.6 0.06],'HorizontalAlignment','left',...
   'BackgroundColor',AQUAMARINE,'Value',0,...
   'FontWeight','demi','FontSize',[12],...
   'Style','radiobutton','String','Non conservative form ','CallBack',...
   'set(CHOICE(1),value,0);set(CHOICE(2),value,1);call_xor(FORMA(1,:),visible,OFF,OFF);call_xor(FORMA(2,:),visible,ON,ON);[dc_sx,dc_dx]=der_co_nom(0);');
%
% Definizione degli estremi dell'intervallo computazionale.
uicontrol('Units','normalized',...
   'Position',[0 0.58 1 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'ForegroundColor','black','String',...
   'in (          ,          ), with boundary conditions  :',...
   'FontAngle','oblique','BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[12]);
OMEGA(1)=uicontrol('Units','normalized',...
   'Position',[0.06 0.58 0.07 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0',...
   'BackgroundColor',AQUAMARINE,'ForegroundColor','blue',...
   'FontWeight','demi','FontSize',[10],'CallBack',...
   'if str2num(get(OMEGA(2),stringa))<=str2num(get(OMEGA(1),stringa)), gest_errori(2); else smista_coef_bc(get(CHOICE(1),value)), end;');

OMEGA(2)=uicontrol('Units','normalized',...
   'Position',[0.16 0.58 0.07 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'CallBack',...
   'if str2num(get(OMEGA(2),stringa))<=str2num(get(OMEGA(1),stringa)), gest_errori(2); else smista_coef_bc(get(CHOICE(1),value)), end;');
%
% Definizione delle condizioni al contorno
uicontrol('Units','normalized',...
   'Position',[0.06 0.52 1 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'ForegroundColor','magenta','String',...
   'Left End                                       Right End     ',...
   'FontAngle','oblique','BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[12]);
%
[dc_sx,dc_dx]=der_co_nom(0);
% Dirichlet a sinistra (default)
BC_SX(1)=uicontrol('Units','normalized',...
   'Position',[0.02 0.46 0.02 0.03],'HorizontalAlignment','left',...
   'BackgroundColor',AQUAMARINE,'Value',1,...
   'FontWeight','demi','FontSize',[12],...
   'Style','radiobutton','String',' ','CallBack',...
   'set(BC_SX(2),value,0);set(BC_SX(3),value,0);set(BC_SX(1),value,1);'); 
uicontrol('Units','normalized',...
   'Position',[0.07 0.45 0.05 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','u  = ','ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
BCVALUE(1,1) =uicontrol('Units','normalized',...
   'Position',[0.15 0.45 0.3 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
%
% Neumann estremo di sinistra
BC_SX(2)=uicontrol('Units','normalized',...
   'Position',[0.02 0.35 0.02 0.03],'HorizontalAlignment','left',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[12],...
   'Style','radiobutton','String',' ','CallBack',...
   'set(BC_SX(1),value,0);set(BC_SX(2),value,1);set(BC_SX(3),value,0);'); 
DER_CON(1,1)=uicontrol('Units','normalized',...
   'Position',[0.05 0.34 0.2 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',cat(2,dc_sx,'= '),'ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
BCVALUE(2,1) =uicontrol('Units','normalized',...
   'Position',[0.25 0.34 0.2  0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
%
% Robin a sinistra
BC_SX(3)=uicontrol('Units','normalized',...
   'Position',[0.02 0.24 0.02 0.03],'HorizontalAlignment','left',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[12],...
   'Style','radiobutton','String',' ','CallBack',...
   'set(BC_SX(2),value,0);set(BC_SX(1),value,0);set(BC_SX(3),value,1);'); 
DER_CON(1,2)=uicontrol('Units','normalized',...
   'Position',[0.05 0.23 0.1 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',cat(2,dc_sx,' + '),'ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
BCVALUE(3,1)=uicontrol('Units','normalized',...
   'Position',[0.15 0.23 0.1 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
uicontrol('Units','normalized',...
   'Position',[0.25 0.23 0.05 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',' u = ','ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
BCVALUE(4,1) =uicontrol('Units','normalized',...
   'Position',[0.3 0.23 0.15 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
%
% Dirichlet a destra (default)
BC_DX(1)=uicontrol('Units','normalized',...
   'Position',[0.52 0.46 0.02 0.03],'HorizontalAlignment','left',...
   'BackgroundColor',AQUAMARINE,'Value',1,...
   'FontWeight','demi','FontSize',[12],...
   'Style','radiobutton','String',' ','CallBack',...
   'set(BC_DX(2),value,0);set(BC_DX(3),value,0);set(BC_DX(1),value,1);'); 
uicontrol('Units','normalized',...
   'Position',[0.57 0.45 0.05 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','u  =      ','ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
BCVALUE(1,2) =uicontrol('Units','normalized',...
   'Position',[0.65 0.45 0.3 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
%
% Neumann a destra
BC_DX(2)=uicontrol('Units','normalized',...
   'Position',[0.52 0.35 0.02 0.03],'HorizontalAlignment','left',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[12],...
   'Style','radiobutton','String',' ','CallBack',...
   'set(BC_DX(1),value,0);set(BC_DX(2),value,1);set(BC_DX(3),value,0);'); 
DER_CON(2,1)=uicontrol('Units','normalized',...
   'Position',[0.55 0.34 0.2 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',cat(2,dc_dx,'= '),'ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
BCVALUE(2,2) =uicontrol('Units','normalized',...
   'Position',[0.75 0.34 0.2 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','0','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
%
% Robin a destra
BC_DX(3)=uicontrol('Units','normalized',...
   'Position',[0.52 0.24 0.02 0.03],'HorizontalAlignment','left',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[12],...
   'Style','radiobutton','String',' ','CallBack',...
   'set(BC_DX(2),value,0);set(BC_DX(1),value,0);set(BC_DX(3),value,1);'); 
DER_CON(2,2)=uicontrol('Units','normalized',...
   'Position',[0.55 0.23 0.1 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',cat(2,dc_dx,' + '),'ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
BCVALUE(3,2)=uicontrol('Units','normalized',...
   'Position',[0.65 0.23 0.1 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
uicontrol('Units','normalized',...
   'Position',[0.75 0.23 0.05 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',' u = ','ForegroundColor','red',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
BCVALUE(4,2) =uicontrol('Units','normalized',...
   'Position',[0.8 0.23 0.15 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
%
% Scelta della modalita' avanzata
MODUS=uicontrol('Units','normalized',...
   'Position',[0.3 0.1 0.425 0.05],'HorizontalAlignment','left',...
   'Style','radiobutton','String','Advanced Mode        ',...
   'FontAngle','oblique','ForegroundColor','red','Value',0,...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[14]);
ESEGUI1=uicontrol('Units','normalized',...
   'Position',[0.1 0.075 0.2 0.1],'HorizontalAlignment','left',...
   'Style','PushButton',...
   'String','GO    ','FontAngle','oblique',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[14],'CallBack',...
   'if str2num(get(OMEGA(2),stringa))<=str2num(get(OMEGA(1),stringa)), gest_errori(2); elseif get(BC_SX(2),value)*get(BC_DX(2),value)==0, femdefell; else, gest_errori(1); end;');
uicontrol('Units','normalized',...
   'Position',[0.75 0.075 0.2 0.1],'HorizontalAlignment','left',...
   'Style','PushButton',...
   'String','Exit','FontAngle','oblique',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[14],'CallBack',...
   'close all; clear all; fem1dmain;');

return

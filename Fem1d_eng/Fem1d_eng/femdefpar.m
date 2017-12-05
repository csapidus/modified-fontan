function femdefpar
%
%FEMDEFELL 1-d finite element definitions.
%
%   See also FEM1D.

%   Reference: A. Quarteroni, Modellistica Numerica 
%   per Problemi Differenziali, Springer-Italia, Milano
%   (2000).

%   Fausto Saleri, Alessandro Veneziani
%   $Revision: 1.1.1.1 $  $Date: 2001/03/09 08:22:48 $
global AQUAMARINE RELDATA METODO SOLVER
global MODUS FEMSYS GERBASE HMESH QUADRA
global FATHER
global THETA DELTAT
%
%[m,n]=size(FATHER)
%if n>2
%   close(FATHER(3));
%end;   
% Posiziona la figura
[figurePos]=pos_figure;
%
% Crea la figura
advanced = get(MODUS,'Value');
if ( advanced == 0) 
   FATHER(3)=figure('NumberTitle','off','Color',AQUAMARINE,...
      'Units','points',...
      'MenuBar','none','Name',...
      'Fem 1D: Finite Element Definition (simple mode).',...
      'Resize','off','Position',figurePos);
else
   FATHER(3)=figure('NumberTitle','off','Color',AQUAMARINE,...
      'Units','points',...
      'MenuBar','none','Name',...
      'Fem 1D: Finite Element Definition (advanced mode).',...
      'Resize','off','Position',figurePos);
end
%
% Costruisce titolo, autori e release
cornici('no ');
uicontrol('Units','normalized',...
   'Position',[0 0.95 1 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','Input the attributes of the numerical scheme ',...
   'BackgroundColor',AQUAMARINE,'FontAngle','oblique',...
   'FontWeight','demi','FontSize',[12]); 
%
% Tipo di elemento finito
uicontrol('Units','normalized',...
   'Position',[0.05 0.8 0.4 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','Finite element type    :',...
   'BackgroundColor',AQUAMARINE,'ForegroundColor','red',...
   'FontWeight','demi','FontSize',[10]); 
s='String';
FEMSYS=uicontrol('Units','normalized',...
   'Position',[0.35 0.81 0.1 0.05],'HorizontalAlignment','left',...
   'Style','Popup','String',['P1';'P2';'P3'],...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10],'ForegroundColor','blue'); 
if advanced ~= 0    % basi gerarchiche solo in modalita' avanzata
   uicontrol('Units','normalized',...
      'Position',[0.47 0.8 0.4 0.05],'HorizontalAlignment','left',...
      'Style','text','String','(hierarchical base         )',...
      'BackgroundColor',AQUAMARINE,'ForegroundColor','red',...
      'FontWeight','demi','FontSize',[10]);
   GERBASE=uicontrol('Units','normalized',...
      'Position',[0.71 0.81 0.02 0.03],'HorizontalAlignment','left',...
      'Style','radiobutton','String',' ',...
      'BackgroundColor',AQUAMARINE,'ForegroundColor','red',...
      'FontWeight','demi','FontSize',[10]);
else
   GERBASE=uicontrol('Units','normalized',...
      'Position',[0.69 0.81 0.02 0.03],'HorizontalAlignment','left',...
      'Style','radiobutton','String',' ','Value',0,...
      'BackgroundColor',AQUAMARINE,'ForegroundColor','red',...
      'FontWeight','demi','FontSize',[10],'Visible','off');
end
% Tipo di griglia
if advanced == 0    % griglie non uniformi solo in modalita' avanzata  
   uicontrol('Units','normalized',...
      'Position',[0.05 0.7 0.4 0.05],'HorizontalAlignment','left',...
      'Style','Text',...
      'String','Grid Spacing  (uniform)     :',...
      'BackgroundColor',AQUAMARINE,'ForegroundColor','red',...
      'FontWeight','demi','FontSize',[10]); 
   HMESH =uicontrol('Units','normalized',...
      'Position',[0.41 0.7 0.1 0.05],'HorizontalAlignment','left',...
      'Style','edit','String','0.1',...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],...
      'ForegroundColor','blue','CallBack','comp_Peclet(0)'); 
   QUADRA =uicontrol('Position',...
      [figurePos(1)+305 figurePos(2)+175 30 20],...
      'Style','edit','String','0',...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],...
      'ForegroundColor','blue','Visible','off');   
else
   uicontrol('Units','normalized',...
      'Position',[0.05 0.7 0.4 0.05],'HorizontalAlignment','left',...
      'Style','Text',...
      'String','Grid Spacing as function of x     :',...
      'BackgroundColor',AQUAMARINE,'ForegroundColor','red',...
      'FontWeight','demi','FontSize',[10]);
   HMESH =uicontrol('Units','normalized',...
      'Position',[0.47 0.7 0.3 0.05],'HorizontalAlignment','left',...
      'Style','edit','String','0.1*(1+exp(-x))',...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','blue');
   uicontrol('Units','normalized',...
      'Position',[0.05 0.6 0.6 0.05],'HorizontalAlignment','left',...
      'Style','Text',...
      'String','N. of nodes in the quadrature formula (GLL):   ',...
      'BackgroundColor',AQUAMARINE,'ForegroundColor','red',...
      'FontWeight','demi','FontSize',[10]);
   QUADRA =uicontrol('Units','normalized',...
      'Position',[0.67 0.6 0.1 0.05],'HorizontalAlignment','left',...
      'Style','edit','String','7',...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','blue');
end
hmesh = get(HMESH,'String');

%
% Tipo di metodo (da scegliersi solo in modalita' avanzata)
if advanced == 0
   METODO=uicontrol('Units','normalized',...
      'Position',[0.05 0.5 0.6 0.05],'HorizontalAlignment','left',...
      'Style','popup','Value',1,...
      'String',['Galerkin              ';...
         'Artificial viscosity  ';...
         'Upwind                ';...
         'Scharfetter-Gummel    ';...
         'Mass-lumping  (FD)    '],...
      'BackgroundColor',AQUAMARINE,'Visible','off',...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','blue');
   SOLVER(1)=uicontrol('Units','normalized',...
      'Position',[0.45 0.41 0.45 0.05],'HorizontalAlignment','left',...
      'Style','popup','Value',1,...
      'String',['LU Factorisation              ';...
         'Cholesky Factorisation        ';...
         'Conjugate Gradient            ';...
         'BiCGStab                      ';...
         'GMRES                         '],...
      'BackgroundColor',AQUAMARINE,'Value',1,'Visible','off',...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','blue'); 
else
   uicontrol('Units','normalized',...
      'Position',[0.05 0.5 0.1 0.05],'HorizontalAlignment','left',...
      'Style','Text',...
      'String','Metodo:',...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','red');
   METODO=uicontrol('Units','normalized',...
      'Position',[0.175 0.5 0.3 0.05],'HorizontalAlignment','left',...
      'Style','popup','Value',1,...
      'String',['Galerkin              ';...
         'Artificial  viscosity ';...
         'Upwind                ';...
         'Scharfetter-Gummel    ';...
         'Mass-lumping (FD)     '],...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','blue'); 
   uicontrol('Units','normalized',...
      'Position',[0.05 0.4 0.4 0.05],'HorizontalAlignment','left',...
      'Style','Text',...
      'String','Risolutore del sistema lineare: ',...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','red');
   SOLVER(1)=uicontrol('Units','normalized',...
      'Position',[0.45 0.41 0.45 0.05],'HorizontalAlignment','left',...
      'Style','popup','Value',1,...
      'String',['LU Factorisation              ';...
         'Conjugate Gradient            ';...
         'BICG                          ';...
         'BiCGStab                      ';...
         'CGS                           ';...
         'GMRES                         '],...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','blue'); 
   uicontrol('Units','normalized',...
      'Position',[0.05 0.3 0.3 0.05],'HorizontalAlignment','left',...
      'Style','Text',...
      'String','Preconditioner    : ',...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','red');
   SOLVER(2)=uicontrol('Units','normalized',...
      'Position',[0.35 0.31 0.45 0.05],'HorizontalAlignment','left',...
      'Style','popup','Value',1,...
      'String',['No  preconditioner ';...
         'Diagonal  norm  2  ';...
         'ILU(0)             ';...
         'MILU               '],...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],'ForegroundColor','blue');
end
uicontrol('Position',[100 30 100 40],'Style','PushButton',...
   'String','GO    ','FontAngle','oblique',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[14],'CallBack',...
   'postprocdefpar');
%
% scelta del parametro theta per il theta metodo
%
uicontrol('Units','normalized',...
   'Position',[0.05 0.2 0.5 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','Time advancement (theta method)    :',...
   'BackgroundColor',AQUAMARINE,'ForegroundColor','black',...
   'FontWeight','demi','FontSize',[10]); 
uicontrol('Units','normalized',...
   'Position',[0.55 0.2 0.4 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','theta = ',...
   'BackgroundColor',AQUAMARINE,'ForegroundColor','red',...
   'FontWeight','demi','FontSize',[10]); 
THETA =uicontrol('Units','normalized',...
      'Position',[0.65 0.2 0.1 0.05],'HorizontalAlignment','left',...
      'Style','edit','String','1',...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],...
      'ForegroundColor','blue','CallBack','th=str2num(get(THETA,stringa));if th>1 | th<0, gest_errori(7);end;'); 
uicontrol('Units','normalized',...
   'Position',[0.55 0.1 0.4 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','Delta t = ',...
   'BackgroundColor',AQUAMARINE,'ForegroundColor','red',...
   'FontWeight','demi','FontSize',[10]); 
DELTAT =uicontrol('Units','normalized',...
      'Position',[0.65 0.1 0.1 0.05],'HorizontalAlignment','left',...
      'Style','edit','String','0.1',...
      'BackgroundColor',AQUAMARINE,...
      'FontWeight','demi','FontSize',[10],...
      'ForegroundColor','blue','CallBack','dt=str2num(get(DELTAT,stringa));if dt<0, gest_errori(9);end;'); 
   
   
return

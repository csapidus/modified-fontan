function pb1d_diffusione_ossigeno
%
% DIFFUSION OF OXYGEN IN A CELLULE
% We study the diffusion of oxygen inside a spherical cellule 
% immersed in an environment with constant concentration.
% The study is made in spherical coordinates, assuming spherical symmetry.
% It is the a one dimensional problem.
%
% Type of problem: Elliptic with constant coefficents (in Cartesian coordinates)
% it becomed with variable coefficents because of the traqnsormation ro spherical 
% coordinates.
%
global AQUAMARINE RELDATA FORMA OMEGA CHOICE BCVALUE
%
% Posiziona la figura
close all;
[figurePos]=pos_figure;
%
% Dichiara il pannello iniziale
father=figure('NumberTitle','off','Color',AQUAMARINE,'Units','points',...
   'MenuBar','none','Name','Fem 1D. Model Problem   ','Resize','off',...
   'Position',figurePos);
%
% Costruisce titolo, autori e release
cornici('no ');
%
% Breve descrizione dei parametri del problema differenziale.
topic={...
   'Consider the diffusion of oxygen into a spherical cellule.'
   '(radius R). Let C_ext be the concentration of oxygen      '
   'at the exterior, that the diffusivity mu is constant and  '
   'that the problem has a spherical symmetry. The law        '
   'governing diffusion is given by                           '
   '                                                          '
   '      - mu Laplace(C) + k C = 0                           '
   '                                                          '
   'Using the spherical symmetry we can formulate the problem '
   'in the radial coordinate. We have:                        '
   '                                                          '
   '  - mu (r^2 C`)` + r^2 k C = 0       0 <= r <= R          '
   '  C`(0)=0; C(R)=C_ext                                     '
   '                                                          '};
uicontrol('Units','normalized',...
   'Position',[0.01 0.65 0.65 0.3],'HorizontalAlignment','left',...
   'Style','list','String',topic,...
   'ForegroundColor','black',... 
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','normal','FontSize',[10]);
%
% Definizione dei parametri del problema modello
uicontrol('Units','normalized',...
   'Position',[0.01 0.575 0.65 0.05],'HorizontalAlignment','left',...
   'Style','Text','ForegroundColor','magenta',...
   'String',...
   'Input problem data                 : ',...
   'BackgroundColor',AQUAMARINE,'FontAngle','oblique',...
   'FontWeight','demi','FontSize',[12]);
%
% Forma dell'equazione
CHOICE(1)=uicontrol('Visible','off','Value',1);
%
% Lunghezza del cuscinetto
uicontrol('Units','normalized',...
   'Position',[0.03 0.475 0.4 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'ForegroundColor','black','String',...
   'Cellule Radius           R = ',...
   'FontAngle','oblique','BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
OMEGA(1)=uicontrol('Visible','off','String','0');
OMEGA(2)=uicontrol('Units','normalized',...
   'Position',[0.45 0.475 0.15 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
% Non capisco perche' a me non funziona!
%,'CallBack',...
%   'if str2num(get(OMEGA(2),stringa))<= 0, gest_errori(2); else smista_coef_bc(get(CHOICE(1),value)), end;');
%
% Diffusivita'
uicontrol('Units','normalized',...
   'Position',[0.03 0.375 0.45 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'ForegroundColor','black','String',...
   'Diffusivity''  mu = ',...
   'FontAngle','oblique','BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
FORMA(1,1)=uicontrol('Units','normalized',...
   'Position',[0.45 0.375 0.3 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1',...
   'BackgroundColor',AQUAMARINE,'ForegroundColor','blue',...
   'FontWeight','demi','FontSize',[10]);
%
FORMA(1,2)=uicontrol('Units','normalized',...
   'Position',[0.45 0.275 0.1 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','0','Visible','off',...
   'BackgroundColor',AQUAMARINE,'ForegroundColor','blue',...
   'FontWeight','demi','FontSize',[10]);
%
uicontrol('Units','normalized',...
   'Position',[0.03 0.275 0.45 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'ForegroundColor','black','String',...
   'Consumption coefficent   k = ',...
   'FontAngle','oblique','BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
FORMA(1,3)=uicontrol('Units','normalized',...
   'Position',[0.45 0.275 0.1 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1',...
   'BackgroundColor',AQUAMARINE,'ForegroundColor','blue',...
   'FontWeight','demi','FontSize',[10]);
FORMA(1,4)=uicontrol('String','0','Visible','off');
uicontrol('Units','normalized',...
   'Position',[0.03 0.175 0.15 0.05],'HorizontalAlignment','left',...
   'Style','Text',...
   'String','C_ext =',...
   'BackgroundColor',AQUAMARINE,'FontAngle','Oblique',...
   'FontWeight','demi','FontSize',[10]);
BCVALUE(1,2) =uicontrol('Units','normalized',...
   'Position',[0.45 0.175 0.3 0.05],'HorizontalAlignment','left',...
   'Style','Edit',...
   'String','1','ForegroundColor','blue',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[10]);
%
% Descrizione grafica del caso test
subplot(3,3,3); theta=[0:0.01:2*pi];
x=cos(theta);y=sin(theta);
%f=[f 1]; x=[x 1];
patch(x,y,'blue'); % patch([0 1 1 0],[-0.3 -0.3 -0.1 -0.1],'red');
testo=text(0,1.2,'Cext');
set(testo,'FontSize',12)
testo=text(0.75,0.2,'R');
set(testo,'FontSize',12);
hold on; 
theta=[0:0.5:2*pi];
x=cos(theta);y=sin(theta);
u=-0.5*cos(theta);v=-0.5*sin(theta);
H=quiver(x,y,u,v,'r');
axis off; hold off;
%
% Crea il bottone di esecuzione
esegui=uicontrol('Units','normalized',...
   'Position',[0.035 0.085 0.25 0.07],'HorizontalAlignment','left',...
   'Style','PushButton',...
   'String','  GO  ','FontAngle','oblique',...
   'BackgroundColor',AQUAMARINE,...
   'FontWeight','demi','FontSize',[14],...
   'CallBack','input2diffox');




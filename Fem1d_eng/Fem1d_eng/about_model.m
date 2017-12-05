function about_model
%
%ABOUT_MODEL 
%   create a box for displaying the 'about' information

%   Fausto Saleri, Alessandro Veneziani
%   $Revision: 1.1.1.1 $  $Date: 2001/03/09 08:22:48 $

global MODELLI MODELS
t=help(MODELS(get(MODELLI,'Value')).name);
uicontrol('Units','normalized',...
   'Position',[0.15 0.2 0.7 0.3],'HorizontalAlignment','left',...
   'Style','Text',...
   'String',t,'ForegroundColor','blue',...
   'BackgroundColor',[127/255 1 1],...
   'FontWeight','demi','FontSize',[8]);
return
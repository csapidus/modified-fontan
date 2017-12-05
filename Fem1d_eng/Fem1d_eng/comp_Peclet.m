function [Pe_l_t,Pe_l_r]=comp_Peclet(advanced)
global AQUAMARINE CHOICE FORMA OMEGA HMESH METODO

a=str2num(get(OMEGA(1),'String'));
b=str2num(get(OMEGA(2),'String'));
if a >= b
   disp('Error at domain boundary'); figure(1);
else
   if get(CHOICE(1),'Value') == 1
      i_forma = 1;
   else
      i_forma = 2;
   end
   diffusion = [get(FORMA(i_forma,1),'String'),'+0.*x'];
   transport = [get(FORMA(i_forma,2),'String'),'+0.*x'];
   reaction = [get(FORMA(i_forma,3),'String'),'+0.*x'];
   h = [get(HMESH,'String'),'+0.*x'];
   x = a; xx(1) = a; k = 1;
   while xx(end) < b
      k = k + 1; xx(k) = xx(k-1) + eval(h);
      x = xx(end);
   end
   xx(end) = b;
   x = xx;
   mu   = eval(diffusion);
   beta = abs(eval(transport));
   sigma = abs(eval(reaction));
   h = eval(h);
   Pe_l_t=0.5*h.*beta./mu;
   Pe_l_r=(h.^2).*sigma./(6*mu);
   Peclet_locale_t = max(Pe_l_t);
   Peclet_locale_r = max(Pe_l_r); stampa = 0;
   stampa = 0;
   if (Peclet_locale_t > Peclet_locale_r)&(Peclet_locale_t > 1)
      Peclet = num2str(fix(Peclet_locale_t));
      stringa1=['COnvection dominated problem, with Peclet = ',Peclet];
      h = 2*mu/beta; h = num2str(h);
      stringa2=['No stabilisation needed if  h < ',h];
      if advanced == 0, set(METODO,'Value',4); stampa = 1;messaggio=5; end      
   elseif Peclet_locale_r > 1
      Peclet = num2str(fix(Peclet_locale_r));
      stringa1=['Reaction dominated problem with Peclet = ',Peclet];
      h = sqrt(6*mu/sigma); h = num2str(h);
      stringa2=['No stabilisation needed if h < ',h];
      if advanced == 0, set(METODO,'Value',5); stampa = 1;messaggio=6; end
   end
   if stampa == 1
      gest_errori(messaggio,h);
%      uicontrol('Units','normalized',...
%         'Position',[0.2 0.25 0.3 0.05],'HorizontalAlignment','left',...
%         'Style','Text','String','Attenzione:',...
%         'BackgroundColor',AQUAMARINE,'FontAngle','oblique',...
%         'FontWeight','demi','FontSize',[10],...
%         'ForegroundColor','magenta');
%      if advanced == 0
%         uicontrol('Units','normalized',...
%            'Position',[0.25 0.2 0.7 0.05],'HorizontalAlignment','left',...
%            'Style','Text','String',stringa1,...
%            'BackgroundColor',AQUAMARINE,...
%            'FontWeight','demi','FontSize',[10],...
%            'ForegroundColor','magenta');
%         uicontrol('Units','normalized',...
%            'Position',[0.25 0.15 0.6 0.05],'HorizontalAlignment','left',...
%            'Style','Text','String',stringa2,...
%            'BackgroundColor',AQUAMARINE,...
%            'FontWeight','demi','FontSize',[10],...
%            'ForegroundColor','magenta');
%      else
%         uicontrol('Units','normalized',...
%            'Position',[0.25 0.2 0.7 0.05],'HorizontalAlignment','left',...
%            'Style','Text','String',stringa1,...
%            'BackgroundColor',AQUAMARINE,...
%            'FontWeight','demi','FontSize',[10],...
%            'ForegroundColor','magenta');
%      end
   end
end

return
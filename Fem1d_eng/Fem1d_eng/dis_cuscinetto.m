function dis_cuscinetto

global FORMA
global uh coord

subplot('Position',[0.05 0.25 0.9 0.42]);
s = get(FORMA(1,1),'String');
[xs,s] = fplot(s,[coord(1),coord(end)]);
plot(coord,uh,'r',xs,s,'--');
legend('Pression','Lubricat. ring surf.');

return
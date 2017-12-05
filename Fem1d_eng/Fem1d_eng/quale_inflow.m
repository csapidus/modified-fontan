function infl=quale_inflow
%
% selezione del bordo di inflow
% infl=0 => inflow a sinistra
% infl=1 => inflow a destra
% infl=2 => inflow da entrambe le parti
% (infl=-1 => non ci sono inflow <<< attualmente non gestito)
%
global AQUAMARINE RELDATA CHOICE METODO
global CHOICE BC_DX BC_SX HMESH QUADRA SOLVER
global MODUS OMEGA FORMA BCVALUE FEMSYS GERBASE DER_CON
global USCITA AN_ERR
global FATHER ON OFF
global MODELLI MODELS SERVICE
%
quale_forma=get(CHOICE(1),'value');
in_sx=0; in_dx=0;
if quale_forma==0, quale_forma=2; end; %non conservativo
x=str2num(get(OMEGA(1),'String'));
a_sx=eval(get(FORMA(quale_forma,1),'String'))
if a_sx > 0, in_sx=1, else in_sx=0, end; %estermo sinistro=inflow
x=str2num(get(OMEGA(2),'String'));
a_dx=eval(get(FORMA(quale_forma,1),'String'))
if a_dx < 0, in_dx=1, else in_dx=0, end; %estermo destro=inflow
in_sx,in_dx
if exist('SERVICE')
if in_sx==1 & in_dx==0
% Definizione delle condizioni al contorno
 set(SERVICE(1),'Visible','on');
 set(SERVICE(2),'Visible','off');
 set(SERVICE(3),'Visible','on');
 %set(SERVICE(4),'Visible','on');
 set(BC_SX(1),'Visible','on');
 set(BC_SX(1),'Value',1);
 set(BCVALUE(1,1),'Visible','on');
% set(BCVALUE(2,1),'Visible','on');
% set(BCVALUE(3,1),'Visible','on');
% set(BCVALUE(4,1),'Visible','on');
 %set(BC_SX(2),'Visible','on');
 %set(DER_CON(1,1),'Visible','on');
 %set(DER_CON(1,2),'Visible','on');
 %set(BC_SX(3),'Visible','on');
 %set(DER_CON(1,2),'Visible','on');
 set(SERVICE(5),'Visible','off');
 %set(SERVICE(6),'Visible','off');
 set(BC_DX(1),'Visible','off');
 set(BC_DX(1),'Value',0);
 set(BCVALUE(1,2),'Visible','off');
 %set(BCVALUE(2,2),'Visible','off');
 %set(BCVALUE(3,2),'Visible','off');
 %set(BCVALUE(4,2),'Visible','off');
 %set(BC_DX(2),'Visible','off');
 %set(DER_CON(2,1),'Visible','off');
 %set(DER_CON(2,2),'Visible','off');
 %set(BC_DX(3),'Visible','off');
 infl=0; 
elseif in_sx==0 & in_dx==1
 set(SERVICE(1),'Visible','off');
 set(SERVICE(2),'Visible','on');
 set(SERVICE(3),'Visible','off');
 %set(SERVICE(4),'Visible','off');
 set(BC_SX(1),'Visible','off');
 set(BC_SX(1),'Value',0);
 set(BCVALUE(1,1),'Visible','off');
 %set(BCVALUE(2,1),'Visible','off');
 %set(BCVALUE(3,1),'Visible','off');
 %set(BCVALUE(4,1),'Visible','off');
 %set(BC_SX(2),'Visible','off');
 %set(DER_CON(1,1),'Visible','off');
 %set(DER_CON(1,2),'Visible','off');
 %set(BC_SX(3),'Visible','off');
 %set(DER_CON(1,2),'Visible','off');
 set(SERVICE(5),'Visible','on');
 %set(SERVICE(6),'Visible','on');
 set(BC_DX(1),'Visible','on');
 set(BC_DX(1),'Value',1);
 set(BCVALUE(1,2),'Visible','on');
 %set(BCVALUE(2,2),'Visible','on');
 %set(BCVALUE(3,2),'Visible','on');
 %set(BCVALUE(4,2),'Visible','on');
 %set(BC_DX(2),'Visible','on');
 %set(DER_CON(2,1),'Visible','on');
 %set(DER_CON(2,2),'Visible','on');
 %set(BC_DX(3),'Visible','on');
 infl=1; 
elseif in_sx==1 & in_dx==1
 set(SERVICE(1),'Visible','on');
 set(SERVICE(2),'Visible','on');
 set(SERVICE(3),'Visible','on');
 %set(SERVICE(4),'Visible','on');
 set(BC_SX(1),'Visible','on');
 set(BC_SX(1),'Value',1);
 set(BCVALUE(1,1),'Visible','on');
 %set(BCVALUE(2,1),'Visible','on');
 %set(BCVALUE(3,1),'Visible','on');
 %set(BCVALUE(4,1),'Visible','on');
 %set(BC_SX(2),'Visible','on');
 %set(DER_CON(1,1),'Visible','on');
 %set(DER_CON(1,2),'Visible','on');
 %set(BC_SX(3),'Visible','on');
 %set(DER_CON(1,2),'Visible','on');
 set(SERVICE(5),'Visible','on');
% set(SERVICE(6),'Visible','on');
 set(BC_DX(1),'Visible','on');
 set(BC_DX(1),'Value',1);
 set(BCVALUE(1,2),'Visible','on');
 %set(BCVALUE(2,2),'Visible','on');
 %set(BCVALUE(3,2),'Visible','on');
 %set(BCVALUE(4,2),'Visible','on');
 %set(BC_DX(2),'Visible','on');
 %set(DER_CON(2,1),'Visible','on');
 %set(DER_CON(2,2),'Visible','on');
 %set(BC_DX(3),'Visible','on');
 infl=2; 
end;

[dc_sx,dc_dx]=der_co_nom(0);
end;
return;

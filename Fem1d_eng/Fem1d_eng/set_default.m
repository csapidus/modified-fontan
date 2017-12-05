function set_default(input)
%
%
%
%
%
  global AQUAMARINE RELDATA CHOICE METODO
  global CHOICE BC_DX BC_SX HMESH QUADRA SOLVER 
  global MODUS OMEGA FORMA BCVALUE FEMSYS GERBASE DER_CON
  global USCITA AN_ERR INFO INITC THETA TIME DELTAT
  global FATHER ON OFF SUPP
  global MODELLI MODELS SERVICE LUMPING
  global coord uh mat_fem coord_old uh_old stringa_uex stringa_uxex
  global infl
% Variabili d'ambiente da non cancellare perche' richiamate in esecuzione
  value='Value'; 
  stringa='String'; 
  visible='Visible'; 
  ON = 'on'; OFF = 'off';
  zero='0'; astra = '****';

  switch input
   case 1 %elliptic 1
    set(FEMSYS,value,1);
  otherwise   
    disp('pippo2');
  end; 
   
  return;
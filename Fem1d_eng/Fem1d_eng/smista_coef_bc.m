function smista_coef_bc(forma,problema)
%
% funzioncina di servizio per dare i coefficienti alle b.c.
% forma = 0 => non conservativa
% forma = 1 => conservativa
%
% problema (opzionale) specifica il tipo di problema
%
global FORMA DER_CON BCVALUE
stringa = 'String';
[dc_sx,dc_dx]=der_co_nom(forma)
if (nargin==2) & (problema~=3) %3= problema iperbolico = solo Dirichlet
set(DER_CON(1,1),stringa,cat(2,dc_sx,'  = '));
set(DER_CON(1,2),stringa,cat(2,dc_sx,'  + '));
set(BCVALUE(2,1),stringa,'0');
set(BCVALUE(3,1),stringa,'1');
set(BCVALUE(4,1),stringa,'0');
set(DER_CON(2,1),stringa,cat(2,dc_dx,'  = '));
set(DER_CON(2,2),stringa,cat(2,dc_dx,'  + '));
set(BCVALUE(2,2),stringa,'0');
set(BCVALUE(3,2),stringa,'1');
set(BCVALUE(4,2),stringa,'0');
end;
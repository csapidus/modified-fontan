function [uhdxx] = stimas2(ne,coord,uh)
for ie = 1:ne
   uhdx(ie) = (uh(ie+1)-uh(ie))/(coord(ie+1)-coord(ie));
   h(ie) = coord(ie+1)-coord(ie);   
end
uhdxx(1)=2*h(1)*(uhdx(2)-uhdx(1))/(h(1)+h(2));
for ie = 2:ne-1
   uhdxx(ie) = h(ie)*(uhdx(ie+1)-uhdx(ie-1))/(h(ie)+h(ie+1)/2+h(ie-1)/2);
end
uhdxx(ne)=2*h(ne)*(uhdx(ne)-uhdx(ne-1))/(h(ne)+h(ne-1));
uhdxx = abs(uhdxx);
return


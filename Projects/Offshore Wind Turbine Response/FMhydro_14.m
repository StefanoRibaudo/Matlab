function [F,M]=FMhydro_14(x0dot,tetadot,t,u,udot)
%passare x0dot,tetadot come numeri e non vettori. MAI passare T=0

rho=1025;
D=9.4;
A=D^2/4*pi;
zbot=-113.4;
Cm=1;
CD=0.6;

z=linspace(zbot,0,114);
xdot=dispx(z,x0dot,tetadot);

%ATTENZIONE: NON SONO COMRPESI I VALORI DIPENDENTI DA XDOTDOT SICCOME SONO
%COMPRESI NELLA MATRICE A
dF=rho*(Cm+1)*A*udot(t,:)+0.5*rho*CD*D*(u(t,:)-xdot).*abs(u(t,:)-xdot);
F=trapz(z,dF);
M=trapz(z,dF.*z);
end
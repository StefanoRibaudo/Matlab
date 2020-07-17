function [F,M]=FMhydro(H,T,x0,teta,x0dot,tetadot,t,phi,k)
%passare x0dot,tetadot come numeri e non vettori. MAI passare T=0
rho=1025;
D=9.4;
A=D^2/4*pi;
zbot=-113.4;
Cm=1;
CD=0.6;

z=linspace(zbot,0,114);
f=1./T;
a=H/2;
w=2*pi*f;
h=320;

x=dispx(z,x0,teta);
u=0;
udot=0;
for i=1:length(a)
    u=u+w(i)*a(i)*cosh(k(i)*(z+h))/sinh(k(i)*h).*cos(w(i)*t-k(i)*x+phi(i));
    udot=udot-w(i)^2*a(i)*cosh(k(i)*(z+h))/sinh(k(i)*h).*sin(w(i)*t-k(i)*x+phi(i));
end
xdot=dispx(z,x0dot,tetadot);


%ATTENZIONE: NON SONO COMRPESI I VALORI DIPENDENTI DA XDOTDOT SICCOME SONO
%COMPRESI NELLA MATRICE A
dF=rho*(Cm+1)*A*udot+0.5*rho*CD*D*(u-xdot).*abs(u-xdot);
F=trapz(z,dF);
M=trapz(z,dF.*z);
end
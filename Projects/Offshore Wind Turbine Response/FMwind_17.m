function [F,M]=FMwind_17(Vhub,x0dot,tetadot,Ct)
rho=1.29; %Kg/m^3
Vrated=11.4; %m/s
Drotor=126; %m

zhub=90; %m
Arotor=0.25*pi*Drotor^2; %m^2
xdot=dispx(zhub,x0dot,tetadot);
Vrel=Vhub-xdot;

F=0.5*rho*Arotor*Ct*Vrel*abs(Vrel);
M=F*zhub;

end
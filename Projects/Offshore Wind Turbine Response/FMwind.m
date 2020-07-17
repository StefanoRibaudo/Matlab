function [F,M]=FMwind(Vhub,x0dot,tetadot)
rho=1.29; %Kg/m^3
Vrated=11.4; %m/s
Drotor=126; %m
Ct0=0.75;
a=0.25;
b=0.86;
zhub=90; %m
Arotor=0.25*pi*Drotor^2; %m^2
xdot=dispx(zhub,x0dot,tetadot);
Vrel=Vhub-xdot;
if Vrel<=Vrated
    Ct=Ct0;
else
    Ct=Ct0*exp(-a*(Vrel-Vrated)^b);
end
F=0.5*rho*Arotor*Ct*Vrel*abs(Vrel);
M=F*zhub;

end
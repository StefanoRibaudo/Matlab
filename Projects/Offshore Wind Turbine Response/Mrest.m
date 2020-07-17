function M=Mrest(x0,teta,zcmtot,mtot)
rho=1025; %Kg/m^3
g=9.81; %m/s^2
D=9.4; %m
zbot=-113.4; %m
zCB=zbot/2;
V=abs(zbot*pi/4*D^2);
FB=rho*g*V;
FG=-mtot*g; %segno meno causa convenzione dei segni, positivo se verso l'alto
MB=-FB*(dispx(zCB,x0,teta)-x0); %segno meno causa: momento positivo se orario
MG=-FG*(dispx(zcmtot,x0,teta)-x0);
M=MB+MG;
end
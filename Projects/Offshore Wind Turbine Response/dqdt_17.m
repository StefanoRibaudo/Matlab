function dqdt = dqdt_17(t,q,M,A,C,B,Vhub,u,udot,tspan2,gamma)
global tindex
t=round(t,2);
for i=tindex:length(tspan2)
    if t==tspan2(i) 
        t=i;
        tindex=i;
        break;
    end
end

[Fh,Mh]=FMhydro_14(q(3),q(4),t,u,udot);
zhub=90; %m
if size(Vhub)==1
    [Fw,Mw]=FMwind_17(Vhub,q(3),q(4),q(5));
    Vrel=Vhub-dispx(zhub,q(3),q(4));
else
    [Fw,Mw]=FMwind_17(Vhub(t),q(3),q(4),q(5));
    Vrel=Vhub(t)-dispx(zhub,q(3),q(4));
end
Vrated=11.4; %m/s
Ct0=0.75;
a=0.25;
b=0.86;
if Vrel<=Vrated
    Ct_the=Ct0;
else
    Ct_the=Ct0*exp(-a*(Vrel-Vrated)^b);
end
dqdt(1)=q(3);
dqdt(2)=q(4);
q3dotq4dot=inv(M+A)*(-C*[q(1);q(2)]-[B*q(3);0]+[Fh+Fw;Mh+Mw]);
dqdt(3)=q3dotq4dot(1);
dqdt(4)=q3dotq4dot(2);
dqdt(5)=-gamma*(q(5)-Ct_the);
dqdt=dqdt';


end
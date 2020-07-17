function dqdt = dqdt_14(t,q,M,A,C,B,Vhub,u,udot,tspan2)
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
if size(Vhub)==1
    [Fw,Mw]=FMwind(Vhub,q(3),q(4));
else
    [Fw,Mw]=FMwind(Vhub(t),q(3),q(4));
end
dqdt(1)=q(3);
dqdt(2)=q(4);
q3dotq4dot=inv(M+A)*(-C*[q(1);q(2)]-[B*q(3);0]+[Fh+Fw;Mh+Mw]);
dqdt(3)=q3dotq4dot(1);
dqdt(4)=q3dotq4dot(2);
dqdt=dqdt';


end
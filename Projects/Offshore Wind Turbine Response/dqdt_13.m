function dqdt = dqdt_13(t,q,M,A,C,B,H,T,Vrel,phi,k)

[Fh,Mh]=FMhydro(H,T,q(1),q(2),q(3),q(4),t,phi,k);
[Fw,Mw]=FMwind(Vrel,q(3),q(4));
dqdt(1)=q(3);
dqdt(2)=q(4);
q3dotq4dot=inv(M+A)*(-C*[q(1);q(2)]-[B*q(3);0]+[Fh+Fw;Mh+Mw]);
dqdt(3)=q3dotq4dot(1);
dqdt(4)=q3dotq4dot(2);
dqdt=dqdt';

end
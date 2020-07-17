function dqdt = dqdt_11(t,q,M,A,C,B,H,T,k)
[F,Mo]=FMhydro(H,T,q(1),q(2),q(3),q(4),t,0,k);


dqdt(1)=q(3);
dqdt(2)=q(4);
q3dotq4dot=inv(M+A)*(-C*[q(1);q(2)]-[B*q(3);0]+[F;Mo]);
dqdt(3)=q3dotq4dot(1);
dqdt(4)=q3dotq4dot(2);
dqdt=dqdt';


end
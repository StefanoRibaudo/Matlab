function dqdt = dqdt_5(t,q,M,A,C,B)

dqdt(1)=q(3);
dqdt(2)=q(4);
q3dotq4dot=inv(M+A)*(-C*[q(1);q(2)]-[B*q(3);0]);
dqdt(3)=q3dotq4dot(1);
dqdt(4)=q3dotq4dot(2);
dqdt=dqdt';

end
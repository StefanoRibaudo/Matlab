Tdur=660; %s
df=1/Tdur; %Hz
tstep=0.1;
fhighcut=0.5; %Hz
% f=df:df:fHighCut;
tspan=0:tstep:Tdur-0.1; %vettore tempo [s]
mfloat=7466330; %kg
mnace=240000; %kg
mtower=249718; %kg
mrotor=110000; %kg
zcmtower=43.4; %m
zcmfloat=-89.92; %m
Icmtower=121700000; %kg*m^2
Icmfloat=4229230000; %kg*m^2
zhub=90; %m
mtot=mfloat+mnace+mtower+mrotor;
zcmtot=(mfloat*zcmfloat+(mrotor+mnace)*zhub+mtower*zcmtower)/mtot;
Itot_0=Icmfloat+mfloat*zcmfloat^2+(zhub)^2*(mrotor+mnace)+Icmtower+zcmtower^2*mtower;
rho=1025; %Kg/m^3
zbot=-113.4; %m
zCB=zbot/2; %m
zmoor=-70; %m
K=41180; %N/m
g=9.81; %m/s^2
Area=9.4^2/4*pi;
Vsub=9.4^2/4*pi*abs(zbot); %m^3, volume immerso
Cm=1;

%Creo le matrici
M=zeros(2,2);
A=zeros(2,2);
C=zeros(2,2);
M(1,1)=mtot;
M(1,2)=mtot*zcmtot;
M(2,1)=mtot*zcmtot;
M(2,2)=Itot_0;
A(1,1)=-rho*Cm*Area*zbot;
A(1,2)=-0.5*rho*Cm*Area*zbot^2;
A(2,1)=-0.5*rho*Cm*Area*zbot^2;
A(2,2)=-1/3*rho*Cm*Area*zbot^3;
C(1,1)=K;
C(1,2)=K*zmoor;
C(2,1)=K*zmoor;
C(2,2)=g*(rho*Vsub*zCB-mtot*zcmtot)+K*zmoor^2;
B=10^5; %N/(m/s)

%Punto 7
freqnat=sqrt(eig(inv(M+A)*C))/(2*pi);
disp(['natural frequency of surge: ', num2str(freqnat(1)), ' Hz'])
disp(['natural frequency of pitch: ', num2str(freqnat(2)), ' Hz'])
%Punto 9
%Guarda funzione dqdt_5

%Punto 10

sol1 = ode4(@(t,q)dqdt_5(t,q,M,A,C,B),tspan,[1;0;0;0]);
sol2 = ode4(@(t,q)dqdt_5(t,q,M,A,C,B),tspan,[0;0.1;0;0]);
sol(:,:,1)=sol1;
sol(:,:,2)=sol2;
% f = (0:(length(tspan)/2))/(Tdur);
figure
for i=1:2
    subplot(2,2,1+(i-1)*2);
    plot(tspan,sol(:,1,i)) %Plot del surge
    title(['Surge, free decay, initial condition: ' num2str(i)])
    ylabel('Surge [m]')
    xlabel('Time [s]')
    grid on
    grid minor
    subplot(2,2,2+(i-1)*2);
    plot(tspan,sol(:,2,i)) %Plot del pitch
    title(['Pitch, free decay, initial condition: ' num2str(i)])
    ylabel('Pitch [rad]')
    xlabel('Time [s]')
    grid on
    grid minor
end
figure
for i=1:2
    dfFFT=1/tspan(end);
    fFFT=(1:length(tspan))*dfFFT-dfFFT;
    a=abs(fft(sol(:,1,i)))/length(tspan);
    pss=2*a.^2/dfFFT;
    pss(1)=0;
    nTrim=round(length(tspan)/2);
    pss(nTrim:end)=0;
    [~,isurge]=max(pss);
    freqsurge(i)=fFFT(isurge);
    
    a=abs(fft(sol(:,2,i)))/length(tspan);
    psp=2*a.^2/dfFFT;
    psp(1)=0;
    psp(nTrim:end)=0;
    [~,ipitch]=max(psp);
    freqpitch(i)=fFFT(ipitch);
    
    subplot(2,2,i);
    plot(fFFT,pss) %Plot del power spectrum, surge
    title(['PS for Surge, free decay, initial condition: ' num2str(i)])
    xlim([0,0.1])
    ylim([min(pss),1.3*max(pss)])
    ylabel('S surge [m^2/Hz]')
    xlabel('Frequency [Hz]')
    line([freqpitch(i) freqpitch(i)], [min(pss) 1.3*max(pss)],'Color','green');
    line([freqsurge(i) freqsurge(i)], [min(pss) 1.3*max(pss)],'Color','red','LineStyle','--');
    legend('Power spectrum',['Pitch frequency = ' num2str(round(freqpitch(i),3))],['Surge frequency = ' num2str(round(freqsurge(i),3))])
    
    subplot(2,2,2+i);
    plot(fFFT,psp) %Plot del power spectrum, pitch
    title(['PS for Pitch, free decay, initial condition: ' num2str(i)])
    xlim([0,0.1])
    ylim([min(psp),1.3*max(psp)])
    line([freqpitch(i) freqpitch(i)], [min(psp) 1.3*max(psp)],'Color','green');
    line([freqsurge(i) freqsurge(i)], [min(psp) 1.3*max(psp)],'Color','red','LineStyle','--');
    legend('Power spectrum',['Pitch frequency = ' num2str(round(freqpitch(i),3))],['Surge frequency = ' num2str(round(freqsurge(i),3))])
    ylabel('S pitch [rad^2/Hz]')
    xlabel('Frequency [Hz]')
end

%Punto 11

sol1 = ode4(@(t,q)dqdt_11(t,q,M,A,C,B,0,1,1),tspan,[1;0;0;0]);
sol2 = ode4(@(t,q)dqdt_11(t,q,M,A,C,B,0,1,1),tspan,[0;0.1;0;0]);
sol(:,:,1)=sol1;
sol(:,:,2)=sol2;
figure
for i=1:2
    subplot(2,2,1+(i-1)*2);
    plot(tspan,sol(:,1,i)) %Plot del surge
    title(['Surge, with hydrodynamic forces, initial condition: ' num2str(i)])
    ylabel('Surge [m]')
    xlabel('Time [s]')
    grid on
    grid minor
    subplot(2,2,2+(i-1)*2);
    plot(tspan,sol(:,2,i)) %Plot del pitch
    title(['Pitch, with hydrodynamic forces, initial condition: ' num2str(i)])
    ylabel('Pitch [rad]')
    xlabel('Time [s]')
    grid on
    grid minor
end

%Punto 12

H=6; %m
T=10; %s
k=ksolve2(1/T,320,g);
sol1 = ode4(@(t,q)dqdt_11(t,q,M,A,C,B,H,T,k),tspan,[0;0;0;0]);

figure
subplot(2,1,1);
plot(tspan,sol1(:,1)) %Plot del surge
title(['Surge, with linear wave'])
ylabel('Surge [m]')
xlabel('Time [s]')
grid on
grid minor
subplot(2,1,2);
plot(tspan,sol1(:,2)) %Plot del pitch
title(['Pitch, with linear wave'])
ylabel('Pitch [rad]')
xlabel('Time [s]')
grid on
grid minor

%Punto 13

sol1 = ode4(@(t,q)dqdt_13(t,q,M,A,C,B,H,T,8,0,k),tspan,[0;0;0;0]);

figure
subplot(2,1,1);
plot(tspan,sol1(:,1)) %Plot del surge
title(['Surge, with steady wind and linear wave'])
ylabel('Surge [m]')
xlabel('Time [s]')
grid on
grid minor
subplot(2,1,2);
plot(tspan,sol1(:,2)) %Plot del pitch
title(['Pitch, with steady wind and linear wave'])
ylabel('Pitch [rad]')
xlabel('Time [s]')
grid on
grid minor

%Punto 14

global tindex
tindex=1;
h=320; %m

%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaJS=3.3;
Hs=6;
Tp=10;

[f,S]=jonswap(Hs,Tp,df,fhighcut,gammaJS);
amp=sqrt(2*S*df);
w=2*pi*f;
phi=2*pi*rand(1,length(f));
k=zeros(1,length(f));
for i=1:length(f)
    k(i)=ksolve2(f(i),320,g);
end
tspan2=0:tstep/2:Tdur-0.1; %vettore tempo [s]
tspan2=round(tspan2,2);
z=linspace(zbot,0,114);
u=zeros(length(tspan2),length(z));
udot=zeros(length(tspan2),length(z));
eta=zeros(1,length(tspan2));
for i=1:length(amp)
    for time=1:length(tspan2)
        eta(time)=eta(time)+amp(i)*cos(w(i)*tspan2(time)+phi(i));
        u(time,:)=u(time,:)+w(i)*amp(i)*cosh(k(i)*(z+h))/sinh(k(i)*h)*cos(w(i)*tspan2(time)+phi(i));
        udot(time,:)=udot(time,:)-w(i)^2*amp(i)*cosh(k(i)*(z+h))/sinh(k(i)*h)*sin(w(i)*tspan2(time)+phi(i));
    end
end

sol1 = ode4(@(t,q)dqdt_14(t,q,M,A,C,B,8,u,udot,tspan2),tspan,[0;0;0;0]);

figure
subplot(2,1,1);
plot(tspan,sol1(:,1)) %Plot del surge
title(['Surge, with irregular waves and steady wind'])
ylabel('Surge [m]')
xlabel('Time [s]')
grid on
grid minor

subplot(2,1,2);
plot(tspan,sol1(:,2)) %Plot del pitch
title(['Pitch, with irregular waves and steady wind'])
ylabel('Pitch [rad]')
xlabel('Time [s]')
grid on
grid minor

%Plot figo con le frequenze

figure

dfFFT=1/tspan(end);
fFFT=(1:length(tspan))*dfFFT-dfFFT;
a=abs(fft(sol1(:,1)))/length(tspan);
pss=2*a.^2/dfFFT;
pss(1)=0;
nTrim=round(length(tspan)/2);
pss(nTrim:end)=0;

a=abs(fft(sol1(:,2)))/length(tspan);
psp=2*a.^2/dfFFT;
psp(1)=0;
psp(nTrim:end)=0;

dfFFT2=1/tspan2(end);
fFFT2=(1:length(tspan2))*dfFFT2-dfFFT2;
a=abs(fft(eta))/length(tspan2);
pse=2*a.^2/dfFFT2;
pse(1)=0;
pse(nTrim:end)=0;
[~,ieta]=max(pse);
freqeta=fFFT2(ieta);

subplot(3,1,1);
plot(fFFT2,pse) %Plot del power spectrum, eta
title(['PS for eta'])
xlim([0,0.15])
linee=line([freqeta freqeta], [min(pse) 1.3*max(pse)],'Color','cyan');
legend([linee],['Eta frequency = ' num2str(round(freqeta,3))])
ylabel('S eta [m^2/Hz]')
xlabel('Frequency [Hz]')

subplot(3,1,2);
plot(fFFT,pss) %Plot del power spectrum, surge
title(['PS for Surge'])
xlim([0,0.15])
lines=line([freqsurge(2) freqsurge(2)], [min(pss) 1.3*max(pss)],'Color','red');
linee=line([freqeta freqeta], [min(pss) 1.3*max(pss)],'Color','cyan');
legend([lines linee],['Surge natural frequency = ' num2str(round(freqsurge(2),3))],['Eta frequency = ' num2str(round(freqeta,3))])
ylabel('S surge [m^2/Hz]')
xlabel('Frequency [Hz]')

subplot(3,1,3);
plot(fFFT,psp) %Plot del power spectrum, pitch
title(['PS for Pitch'])
xlim([0,0.15])
ylabel('S pitch [rad^2/Hz]')
xlabel('Frequency [Hz]')
linep=line([freqpitch(2) freqpitch(2)], [min(psp) 1.3*max(psp)],'Color','green');
linee=line([freqeta freqeta], [min(psp) 1.3*max(psp)],'Color','cyan');
lines=line([freqsurge(2) freqsurge(2)], [min(psp) 1.3*max(psp)],'Color','red');
legend([linep linee lines],['Pitch natural frequency = ' num2str(round(freqpitch(2),3))],['Eta frequency = ' num2str(round(freqeta,3))],['Surge natural frequency = ' num2str(round(freqsurge(2),3))])

%Punto 15
tindex=1;

V10=8;
I=0.14;
l=340.2;
[Vhub,Swind]=Kaimal(V10,I,l,f,tspan2);

sol1 = ode4(@(t,q)dqdt_14(t,q,M,A,C,B,Vhub,u,udot,tspan2),tspan,[0;0;0;0]);

figure
subplot(2,1,1);
plot(tspan,sol1(:,1)) %Plot del surge
title(['Surge, with irregular waves and unsteady wind'])
xlabel('Time [s]')
ylabel('Surge [m]')
grid on
grid minor

subplot(2,1,2);
plot(tspan,sol1(:,2)) %Plot del pitch
title(['Pitch, with irregular waves and unsteady wind'])
xlabel('Time [s]')
ylabel('Pitch [rad]')
grid on
grid minor

%Plot figo con le frequenze

figure

dfFFT=1/tspan(end);
fFFT=(1:length(tspan))*dfFFT-dfFFT;
a=abs(fft(sol1(:,1)))/length(tspan);
pss=2*a.^2/dfFFT;
pss(1)=0;
nTrim=round(length(tspan)/2);
pss(nTrim:end)=0;
[~,isurge]=max(pss);

a=abs(fft(sol1(:,2)))/length(tspan);
psp=2*a.^2/dfFFT;
psp(1)=0;
psp(nTrim:end)=0;
[~,ipitch]=max(psp);

dfFFT2=1/tspan2(end);
fFFT2=(1:length(tspan2))*dfFFT2-dfFFT2;
a=abs(fft(Vhub))/length(tspan2);
psw=2*a.^2/dfFFT2;
psw(1)=0;
nTrim=round(length(tspan2)/2);
psw(nTrim:end)=0;
[~,iwind]=max(psw);
freqwind=fFFT2(iwind);

a=abs(fft(eta))/length(tspan2);
pse=2*a.^2/dfFFT2;
pse(1)=0;
pse(nTrim:end)=0;
[~,ieta]=max(pse);
freqeta=fFFT2(ieta);

subplot(4,1,1);
plot(fFFT2,psw) %Plot del power spectrum, wind speed
title(['PS for Wind speed'])
xlim([0,0.15])
linew=line([freqwind freqwind], [min(psw) 1.3*max(psw)],'Color','magenta');
legend([linew],['Wind frequency = ' num2str(round(freqwind,3))])
ylabel('S wind speed [(m/s)^2/Hz]')
xlabel('Frequency [Hz]')

subplot(4,1,2);
plot(fFFT2,pse) %Plot del power spectrum, eta
title(['PS for eta'])
xlim([0,0.15])
linee=line([freqeta freqeta], [min(pse) 1.3*max(pse)],'Color','cyan');
legend([linee],['Eta frequency = ' num2str(round(freqeta,3))])
ylabel('S eta [m^2/Hz]')
xlabel('Frequency [Hz]')

subplot(4,1,3);
plot(fFFT,pss) %Plot del power spectrum, surge
title(['PS for Surge'])
xlim([0,0.15])
lines=line([freqsurge(2) freqsurge(2)], [min(pss) 1.3*max(pss)],'Color','red');
linee=line([freqeta freqeta], [min(pss) 1.3*max(pss)],'Color','cyan');
linew=line([freqwind freqwind], [min(pss) 1.3*max(pss)],'Color','magenta');
legend([lines linee linew],['Surge natural frequency = ' num2str(round(freqsurge(2),3))],['Eta frequency = ' num2str(round(freqeta,3))],['Wind frequency = ' num2str(round(freqwind,3))])
ylabel('S surge [m^2/Hz]')
xlabel('Frequency [Hz]')

subplot(4,1,4);
plot(fFFT,psp) %Plot del power spectrum, pitch
title(['PS for Pitch'])
xlim([0,0.15])
ylabel('S pitch [rad^2/Hz]')
xlabel('Frequency [Hz]')
linep=line([freqpitch(2) freqpitch(2)], [min(psp) 1.3*max(psp)],'Color','green');
linee=line([freqeta freqeta], [min(psp) 1.3*max(psp)],'Color','cyan');
linew=line([freqwind freqwind], [min(psp) 1.3*max(psp)],'Color','magenta');
legend([linep linee linew],['Pitch natural frequency = ' num2str(round(freqpitch(2),3))],['Eta frequency = ' num2str(round(freqeta,3))],['Wind frequency = ' num2str(round(freqwind,3))])

%Punto 16
H=0;
Tdur=900; %s
df=1/Tdur; %Hz
tspan=0:0.1:Tdur-0.1; %vettore tempo [s]

sol1 = ode4(@(t,q)dqdt_13(t,q,M,A,C,B,H,T,10,0,k),tspan,[0;0;0;0]);
sol2 = ode4(@(t,q)dqdt_13(t,q,M,A,C,B,H,T,16,0,k),tspan,[0;0;0;0]);
figure
subplot(2,2,1);
plot(tspan,sol1(:,1)) %Plot del surge
title(['Surge, V10=10[m/s]'])
subplot(2,2,2);
plot(tspan,sol1(:,2)) %Plot del pitch
title(['Pitch, V10=10[m/s]'])
subplot(2,2,3);
plot(tspan,sol2(:,1)) %Plot del surge
title(['Surge, V10=16[m/s]'])
subplot(2,2,4);
plot(tspan,sol2(:,2)) %Plot del pitch
title(['Pitch, V10=16[m/s]'])

%Punto 17
%Parte 1: condizioni uguali al punto 15
Tdur=660; %s
tspan=0:0.1:Tdur-0.1; %vettore tempo [s]
tindex=1;
gamma=1;
sol1 = ode4(@(t,q)dqdt_17(t,q,M,A,C,B,Vhub,u,udot,tspan2,gamma),tspan,[0;0;0;0;0.75]);
figure
subplot(3,1,1);
plot(tspan,sol1(:,1)) %Plot del surge
title(['Surge, dynamic control, wind-wave climate as Q15'])
xlabel('Time [s]')
ylabel('Surge [m]')
grid on
grid minor

subplot(3,1,2);
plot(tspan,sol1(:,2)) %Plot del pitch
title(['Pitch, dynamic control, wind-wave climate as Q15'])
xlabel('Time [s]')
ylabel('Pitch [rad]')
grid on
grid minor

subplot(3,1,3);
plot(tspan,sol1(:,5)) %Plot di Ct
title(['Ct, dynamic control, wind-wave climate as Q15'])
xlabel('Time [s]')
ylabel('Ct [-]')
ylim([0.5,0.8])
grid on
grid minor

%Parte 2: condizioni come punto 16 b)
Vhub=16;

u=zeros(length(tspan2),length(z));
udot=zeros(length(tspan2),length(z));

tindex=1;
gamma=1;
sol1 = ode4(@(t,q)dqdt_17(t,q,M,A,C,B,Vhub,u,udot,tspan2,gamma),tspan,[0;0;0;0;0.75]);
figure
subplot(3,1,1);
plot(tspan,sol1(:,1)) %Plot del surge
title(['Surge, dynamic control, wind-wave climate as Q16-b)'])
xlabel('Time [s]')
ylabel('Surge [m]')
grid on
grid minor

subplot(3,1,2);
plot(tspan,sol1(:,2)) %Plot del pitch
title(['Pitch, dynamic control, wind-wave climate as Q16-b)'])
xlabel('Time [s]')
ylabel('Pitch [rad]')
grid on
grid minor

subplot(3,1,3);
plot(tspan,sol1(:,5)) %Plot di Ct
title(['Ct, dynamic control, wind-wave climate as Q16-b)'])
xlabel('Time [s]')
ylabel('Ct [-]')
grid on
grid minor

%Punto 18
Vhub=16;
u=zeros(length(tspan2),length(z));
udot=zeros(length(tspan2),length(z));

% gamma=[2 1 0.5 0.3 0.1 0.0001];
gamma=linspace(0.1,1,10);
for i=1:length(gamma)
tindex=1;

sol1 = ode4(@(t,q)dqdt_17(t,q,M,A,C,B,Vhub,u,udot,tspan2,gamma(i)),tspan,[0;0;0;0;0.75]);
figure
subplot(3,1,1);
plot(tspan,sol1(:,1)) %Plot del surge
title(['Surge, dynamic control, wind-wave climate as Q16-b), gamma=' num2str(gamma(i))])
xlabel('Time [s]')
ylabel('Surge [m]')
grid on
grid minor

subplot(3,1,2);
plot(tspan,sol1(:,2)) %Plot del pitch
title(['Pitch, dynamic control, wind-wave climate as Q16-b), gamma=' num2str(gamma(i))])
xlabel('Time [s]')
ylabel('Pitch [rad]')
grid on
grid minor

subplot(3,1,3);
plot(tspan,sol1(:,5)) %Plot di Ct
title(['Ct, dynamic control, wind-wave climate as Q16-b), gamma=' num2str(gamma(i))])
xlabel('Time [s]')
ylabel('Ct [-]')
grid on
grid minor
end
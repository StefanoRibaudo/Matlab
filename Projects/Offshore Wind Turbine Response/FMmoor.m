function [F,M]=FMmoor(x0,teta)
%calcolato assumendo che il corpo ruota attorno al pelo dell'acqua
K=41180; %N/m
zmoor=-70; %m
F=-K*dispx(zmoor,x0,teta);
M=-F*zmoor;

end
function [A,k]=waspmethod(wsp)
% Outputs parameters A and k for the Weibul distribution given wind speeds
% in wsp. This method is based on the third non-central moment and median
% of the sample distribution.

%INPUT:
%wsp = 1-D array of wind speeds. Nans allowed

%OUTPUT:
% A,k = paramenters of Weibul distribution.

options=optimset('Display','off','MaxFunEvals',3000,'TolFun',10^(-10),'MaxIter',3000);
nbins=50;
maxwsp=ceil(max(wsp));
binsdels=linspace(0,maxwsp,nbins+1);
xpdfs=binsdels(1:end-1)+binsdels(2)/2;
pdfs=zeros(1,nbins);
for i=1:length(wsp)
    if ~isnan(wsp(i))
        ind=wsp(i)==binsdels;
        if sum(ind)~=1
            ind=(wsp(i)-maxwsp/nbins)<binsdels & wsp(i)>=binsdels;
        end
        ind(end)=[];
        pdfs(ind)=pdfs(ind)+1;
    end
end
pdfs=pdfs/sum(~isnan(wsp))/binsdels(2);

m3=trapz([0,xpdfs],([0,xpdfs].^3).*[0,pdfs]); %#ok<*NASGU>
P50=median(wsp,'omitnan');
cdfp=0.5;

fun=@(x) [cdfp-exp(-(P50/x(1))^x(2)),x(1)^3*gamma(1+3/x(2))-m3];
sol=fsolve(fun,[1,1],options);
k=sol(2);
A=sol(1);


end









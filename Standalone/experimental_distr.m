function [distr,central_bin_value]=experimental_distr(data,nbins,minvalue,maxvalue)
%INPUT:
%data = any array of input data
%nbins = number of bins for the experimental distribution
%min / maxvalue = min and max values for the distribution
%OUTPUT:
%distr = pdf values (y axis)
%central_bin_value = associated sample space values (x axis)

%example:
%r=normrnd(0,1,[100,100000]);
%[distr,central_bin_value]=experimental_distr(r,100,-3,+3);
%figure
%bar(central_bin_value,distr)

data=data(~isnan(data(:)));
half_binwidth=(maxvalue-minvalue)/(2*nbins);
edges=linspace(minvalue,maxvalue,nbins+1);
central_bin_value=edges(1:end-1)+half_binwidth;
ndata=numel(data);
distr=nan(1,nbins);

for i=1:(nbins-1)
    distr(i)=sum(data(:)>=edges(i) & data(:)<edges(i+1))/ndata/(2*half_binwidth);
end
distr(nbins)=sum(data(:)>=edges(nbins) & data(:)<=edges(nbins+1))/ndata/(2*half_binwidth);

end
% The objective is to compute and plot wind roses for mean wind speed,
% turbulence intensity and probability of occurrence for the data in
% wind_stats.csv. The roses will have 12 sector as per standard.


clearvars
stat=csvread("wind_stats.csv",1,0);
wsp=stat(:,2);
sd=stat(:,3);
dir=stat(:,4);
clearvars stat
fieldnames=["a","b","c","d","e","f","g","h","i","l","m","n"];
for i=1:length(fieldnames)
    sec.wsp.(fieldnames(i))=[];
    sec.sd.(fieldnames(i))=[];
    sec.dir.(fieldnames(i))=[];
end
for i=1:length(dir)
    if dir(i)>345
        dir(i)=dir(i)-360;
    end
    for j=0:30:330
        if dir(i)<j+15 && dir(i)>j-15
            name=fieldnames((j+30)/30);
            sec.wsp.(name)(length(sec.wsp.(name))+1)=wsp(i);
            sec.sd.(name)(length(sec.sd.(name))+1)=sd(i);
            sec.dir.(name)(length(sec.dir.(name))+1)=dir(i);
        end
    end
end
for i=1:length(fieldnames)
    sec.U.(fieldnames(i))=mean(sec.wsp.(fieldnames(i)));
    sec.TI.(fieldnames(i))=sec.sd.(fieldnames(i))./sec.wsp.(fieldnames(i));
    sec.prob.(fieldnames(i))=length(sec.wsp.(fieldnames(i)))/length(wsp);
    binedges(i)=deg2rad(-15+(i-1)*30);
    tempwsp(i)=mean(sec.wsp.(fieldnames(i)));
    tempsd(i)=mean(sec.sd.(fieldnames(i)));
    tempprob(i)=sec.prob.(fieldnames(i));
    
end
binedges(13)=deg2rad(345);
figure
polarhistogram('BinEdges',binedges,'BinCounts',tempwsp)
title("Mean wind speed [m/s]")
rotate_axes()
figure
polarhistogram('BinEdges',binedges,'BinCounts',tempsd)
title("Turbulence intensity")
rotate_axes()
figure
polarhistogram('BinEdges',binedges,'BinCounts',tempprob)
title("Probability of occurrence")
rotate_axes()

function rotate_axes()
pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
end

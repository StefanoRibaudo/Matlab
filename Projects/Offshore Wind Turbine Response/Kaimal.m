function [V,Swind]=Kaimal(V10,I,l,f,time)           % Kaimal wind spectrum function
    %
    Swind=(4*I^2*V10*l)./((1+6*f*l/V10).^(5/3));    % Swind = Kaimal wind spectrum
    %
    b=sqrt(2*Swind*(f(2)-f(1)));                    % b = Wind speed fluctuation amplitude
    %
    e=rand(length(f))*2*pi; %proposta
    for t=1:length(time)                            % Wind time series generator
        Vfluct=0;
        for p=1:length(f)
            %e=rand()*2*pi;
            Vfluct=Vfluct+b(p)*cos(2*pi*f(p)*time(t)+e(p));
        end
        V(t)=V10+Vfluct;
    end

end
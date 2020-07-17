function [f,S]=jonswap(Hs,Tp,df,fHighCut,gammaJS)       % Function for JONSWAP Spectrum computation

    f=df:df:fHighCut; 
    fp=1/Tp;                                            % fp = significant wave frequency for 50-year sea state [Hz];
    
    S=0.3125*Hs^2*Tp*((f/fp).^(-5)).*exp(-1.25*(f/fp).^(-4))*(1-0.287*log(gammaJS)).*gammaJS.^(exp(-0.5*((f/fp-1)/sigma(f,fp)).^2));
                                                        % Spectrum at fp frequency

end
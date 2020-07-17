function k = ksolve2(f,h,g)     % Function for irregular value of k

    kvett=0.0001:0.0001:1;
    eq=(2*pi*f)^2 - g*kvett(1)*tanh(kvett(1)*h);

    flag=0;
    if eq>0
        flagsign=1;
    else
        flagsign=0;
    end
    for i=2:length(kvett)
        eq=(2*pi*f)^2 - g*kvett(i)*tanh(kvett(i)*h);
        if flagsign==1
            if eq<=0
                endi=i;
                flag=1;
                break
            end
        end
        if flagsign==0
            if eq>0
                endi=i;
                flag=1;
                break
            end
        end     
    end

    if flag==1
        k=kvett(endi); 
    else
        k=0.0355;
    end

end


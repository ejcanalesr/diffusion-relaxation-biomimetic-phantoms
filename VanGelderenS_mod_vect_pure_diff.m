function Svg = VanGelderenS_mod_vect_pure_diff(x, bval)
        
        D0     = 2.1; %um2/ms
        Delta  = 35;
        delta  = 9;
        
        %gyroMagnRatio =  267.513 * 10^(-6);
        gyroMagnRatio =  267.5153151 * 10^(-6);

        %bval   =  gyroMagnRatio^2 * G1^2 * delta1^2 * (Delta1 - delta1/3) %ms/um^2
        
        g   = sqrt(bval./( gyroMagnRatio^2 * delta.^2 .* (Delta - delta./3) ));
        q   = g*gyroMagnRatio;

        % q    = sqrt(bval./(Delta - delta/3))./delta;
        
        r    = x(1);
        beta = x(2);
        
        td       = r^2/D0;
        bardelta = delta/td;
        barDelta = Delta/td;

        N=15; 
        b = [1.8412    5.3314    8.5363   11.7060   14.8636   18.0155   21.1644 24.3113   27.4571   30.6019 ...
             33.7462   36.8900   40.0334   43.1766   46.3196 49.4624   52.6050   55.7476   58.8900   62.0323];
   
        s = 0;
        for k=1:N
           s = s + (2/(b(k)^6*(b(k)^2-1)))*(-2 + 2*b(k)^2*bardelta + ...
                          2*(exp(-b(k)^2*bardelta)+exp(-b(k)^2*barDelta)) - ...
                          exp(-b(k)^2*(bardelta+barDelta)) - exp(-b(k)^2*(barDelta-bardelta)));          
        end
        D_b = s.*D0.*(q.^2).*td^3;
        
        D   = D_b./bval;
        Svg = beta .* exp(-D_b) .* sqrt(pi)/2 .* erf(sqrt(bval.*(D0-D)))./sqrt(bval.*(D0-D));
        %Svg = beta * exp(-D_b) .* erf(sqrt(bval.*(D0-D)))./sqrt(bval.*(D0-D));
        %Svg = beta * exp(-D_b) .* 1./sqrt(bval.*(D0-D));
        %Svg = beta * exp(-D_b) .* 1./sqrt(bval);
end
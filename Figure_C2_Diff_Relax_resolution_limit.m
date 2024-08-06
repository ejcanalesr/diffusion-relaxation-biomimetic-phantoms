clear all; clc;

radius = [0.0:0.001:10]; % um

D0     = 2.1; %um2/ms
Delta  = 35;
delta  = 9;
gyroMagnRatio =  267.5153151*10^(-6);
g      = 235.85;
%b_vec =  gyroMagnRatio^2 * g.^2 * delta^2 * (Delta - delta/3); %ms/um^2
b_vec = 10;
q_vec = sqrt(b_vec./(Delta - delta/3))./delta;

T2_CSF = 3000;
rho    = 0.0037;
TE1    = 100; % ms

for i=1:length(radius)
    r = radius(i);
    
    signal_d   = VanGelderenS_mod_vect_pure_diff([r, 1], b_vec); % full signal
    
    signal_r  = VanGelderenS_mod_vect_pure_rel([r, 1], [TE1, T2_CSF, rho]);

    Filter_D(i)  = signal_d;
    Filter_R(i)  = signal_r;    

end

Filter_D     =  Filter_D/max(Filter_D);
Filter_R     =  Filter_R/max(Filter_R);


% ---------------- Plot figure -------------------------------------------%
fig = figure('Renderer', 'painters', 'Position', [100 100 1000 700]);
left_color  = [0 0.4470 0.7410];
right_color = [0.6350 0.0780 0.1840];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

hold on

plot(radius, Filter_D,...
    'color', [0 0.4470 0.7410],...
    'LineWidth',3,...
    'MarkerSize',15,...
    'DisplayName', '\it S(r): van Gelderen, Eqs.(4),(A1) '); 

plot(radius, Filter_R,...
    'color', [0.6350 0.0780 0.1840],...
    'LineWidth',2,...
    'MarkerSize',15,...
    'DisplayName', '\it S(r): T2 relaxation model, Eqs. (2),(3)');

set(gca,'XMinorTick','on','YMinorTick','on');
grid on;
xlabel(' Radius \it(r, \mum)')
ylabel('Spherical mean signal, normalised intensity (a.u.)')

TextFontSize   = 20;
LegendFontSize = 18;

set(0,'DefaultAxesFontName','Times',...
    'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8);
set(gca,'FontName','Times New Roman','FontSize',TextFontSize);

set(gca, 'Box', 'on');

SNR = 100;
thr = 1/SNR;

[val1, ind1] = min(abs(Filter_D - (1-thr)));
Rc1          = radius(ind1);

[val2, ind2] = min(abs(Filter_R - (thr)));
Rc2          = radius(ind2);

plot([1 1] * Rc1, [0 Filter_D(ind1)],   'color', [0.5 0.5 0.5], 'LineStyle', '-.', 'DisplayName', '\it  Diffusion resolution limit (SNR=100)')

hold on

plot([1 1] * Rc2, [0 1], '-.', 'color', [0.4660 0.6740 0.1880], 'DisplayName', ['\it  T2 resolution limit (SNR = 100, TE=100 ms)'])

set(gca, 'XTick', [1:10])

hl = legend('location', 'best');
set(hl, 'interpreter', 'tex')

print(fig, 'Figures/FigureC2_Resolution_limits','-r600','-dpng');

%------------------- Private function ------------------------------------%


function Svg = VanGelderenS_mod_vect_pure_diff(x, bval)
        % A version from the Veraart's function (in Github)
        
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
end

clear all; clc;

radius = [0.1:0.1:10]; % um

D0     = 2.1; %um2/ms
Delta  = 35;
delta  = 9;
gyroMagnRatio =  267.5153151*10^(-6);
g      = 235.85;
%b_vec =  gyroMagnRatio^2 * g.^2 * delta^2 * (Delta - delta/3); %ms/um^2
b_vec = 10;
q_vec = sqrt(b_vec./(Delta - delta/3))./delta;


for i=1:length(radius)
    r = radius(i);
    
    signal_r   = VanGelderenS_mod_vect_pure_diff([r, 1], b_vec); % full signal
    signal_C   = VanGelderen_my_approx_compact(delta, Delta, r, D0, b_vec);
    signal_N   = VanGelderen_Neuman(delta, Delta, r, D0, b_vec);
    signal_NT  = VanGelderen_Neuman_Taylor(delta, Delta, r, D0, b_vec);

    Filter_F(i)  = signal_r;    
    Filter_C(i)  = signal_C;
    Filter_N(i)  = signal_N;
    Filter_NT(i) = signal_NT;

end

Filter_F     =  Filter_F/Filter_F(1);
Filter_C     =  Filter_C/Filter_C(1);
Filter_N     =  Filter_N/Filter_N(1);
Filter_NT    =  Filter_NT/Filter_NT(1);

% ---------------- Plot figure -------------------------------------------%
fig = figure('Renderer', 'painters', 'Position', [100 100 1000 700]);
left_color  = [0 0.4470 0.7410];
right_color = [0.6350 0.0780 0.1840];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

hold on

plot(radius, Filter_F,...
    'color', [0 0.4470 0.7410],...
    'LineWidth',3,...
    'MarkerSize',15,...
    'DisplayName', '\it S(r): van Gelderen, Eqs.(4),(A1) '); 


plot(radius, Filter_C,...
    'color', [0.9290 0.6940 0.1250],...
    'LineWidth',2,...
    'LineStyle', '--',...
    'MarkerSize',15,...
    'DisplayName', '\it S(r): medium-pulse, Eqs. (4),(A3)');

plot(radius, Filter_N,...
    'color', [0.6350 0.0780 0.1840],...
    'LineWidth',2,...
    'MarkerSize',15,...
    'DisplayName', '\it S(r): Neuman long-pulse limit, Eqs.(4),(A2)'); 

ind_val = Filter_NT >= 0;
plot(radius(ind_val), Filter_NT(ind_val),...
    'color', [0 0.6470 0.5410],...
    'LineWidth',2,...
    'MarkerSize',15,...
    'DisplayName', '\it S(r): Neuman-Taylor, exp(-x)=1-x'); 

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

[val_n, ind_n] = min( abs(abs(Filter_F - Filter_N) - 1/100) );
Rn = radius(ind_n);
plot([1 1] * Rn, [0 Filter_F(ind_n)],   'color', [0.5 0.5 0.5], 'LineStyle', '-.', 'DisplayName', '\it  Neuman’s limit of validity')

set(gca, 'XTick', [1:15])

hl = legend('location', 'best');
set(hl, 'interpreter', 'tex')

print(fig, 'Figures/FigureC1_Signal_approximations','-r600','-dpng');

%------------------- Private function ------------------------------------%

function [Svg, D] = VanGelderen_my_approx_compact(delta, Delta, r, D0, bval)

        a = 1.8412;
        
        D = 7/48 * r^4 * (1  - (12/41) * r^2/(D0*delta) * ( 1  - exp(-a^2*D0*delta/r^2) ) ) / ( D0 * delta * (Delta - delta/3) );

        Svg = exp(-D*bval) * erf(sqrt(bval*(D0 - D))) * sqrt(pi/(4* bval * (D0 - D)));
end

function [Svg, D] = VanGelderen_Neuman(delta, Delta, r, D0, bval)
        
        D = 7/48 * r^4 / (D0 * delta * (Delta - delta/3));

        Svg = exp(-D*bval) * sqrt(pi/(4* bval * D0));
end

function [Svg, D] = VanGelderen_Neuman_Taylor(delta, Delta, r, D0, bval)
        
        D = 7/48 * r^4 / (D0 * delta * (Delta - delta/3));

        Svg = (1-D*bval) * sqrt(pi/(4* bval * D0));
end


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

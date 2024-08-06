close all; clear; clc; warning off

% Load GT radius
Phantom1 = load('Histological_radius/ROI1.mat');
Phantom3 = load('Histological_radius/ROI3.mat');
Phantom4 = load('Histological_radius/ROI4.mat');
Phantom5 = load('Histological_radius/ROI5.mat');

Threshold = 40; % discard huge pores and convert to radius

R11 = Phantom1.ROI1.s1; R11 = R11(R11<Threshold)/2;
R12 = Phantom1.ROI1.s2; R12 = R12(R12<Threshold)/2;
R13 = Phantom1.ROI1.s3; R13 = R13(R13<Threshold)/2;
R14 = Phantom1.ROI1.s4; R14 = R14(R14<Threshold)/2;
R15 = Phantom1.ROI1.s5; R15 = R15(R15<Threshold)/2;

R31 = Phantom3.ROI3.s1; R31 = R31(R31<Threshold)/2;
R32 = Phantom3.ROI3.s2; R32 = R32(R32<Threshold)/2;
R33 = Phantom3.ROI3.s3; R33 = R33(R33<Threshold)/2;
R34 = Phantom3.ROI3.s4; R34 = R34(R34<Threshold)/2;
R35 = Phantom3.ROI3.s5; R35 = R35(R35<Threshold)/2;

R41 = Phantom4.ROI4.s1; R41 = R41(R41<Threshold)/2;
R42 = Phantom4.ROI4.s2; R42 = R42(R42<Threshold)/2;
R43 = Phantom4.ROI4.s3; R43 = R43(R43<Threshold)/2;
R44 = Phantom4.ROI4.s4; R44 = R44(R44<Threshold)/2;
R45 = Phantom4.ROI4.s5; R45 = R45(R45<Threshold)/2;

R51 = Phantom5.ROI5.s1; R51 = R51(R51<Threshold)/2;
R52 = Phantom5.ROI5.s2; R52 = R52(R52<Threshold)/2;
R53 = Phantom5.ROI5.s3; R53 = R53(R53<Threshold)/2;
R54 = Phantom5.ROI5.s4; R54 = R54(R54<Threshold)/2;
R55 = Phantom5.ROI5.s5; R55 = R55(R55<Threshold)/2;

R{1} = [R11; R12; R13; R14; R15];
R{2} = [R11; R12; R13; R14; R15];
R{3} = [R31; R32; R33; R34; R35];
R{4} = [R41; R42; R43; R44; R45];
R{5} = [R51; R52; R53; R54; R55];


% Load Spherical mean of dMRI data and standard deviation
SM    = load('SMT_diffusion_all_bs/SMT_all_bs.txt');
SM_dMRI_P1 = SM(6,:); % Matching ROIs numbers/order in images and masks
SM_dMRI_P2 = SM(1,:); 
SM_dMRI_P3 = SM(2,:); 
SM_dMRI_P4 = SM(3,:); 
SM_dMRI_P5 = SM(4,:);

STD_dMRI  = load('SMT_diffusion_all_bs/SMT_std_log_all_bs.txt');
STD_dMRI_P1 = STD_dMRI(6,:);
STD_dMRI_P2 = STD_dMRI(1,:); 
STD_dMRI_P3 = STD_dMRI(2,:); 
STD_dMRI_P4 = STD_dMRI(3,:); 
STD_dMRI_P5 = STD_dMRI(4,:);

SM_S0    = load('SMT_diffusion_all_bs/SMT_all_b0s.txt');
SM_S0_P1 = mean(SM_S0(6,:)); std1 = std(SM_S0(6,:));
SM_S0_P2 = mean(SM_S0(1,:)); std2 = std(SM_S0(1,:));
SM_S0_P3 = mean(SM_S0(2,:)); std3 = std(SM_S0(2,:));
SM_S0_P4 = mean(SM_S0(3,:)); std4 = std(SM_S0(3,:));
SM_S0_P5 = mean(SM_S0(4,:)); std5 = std(SM_S0(4,:));

% Create array structure
SM_P(:,1) = SM_dMRI_P1;
SM_P(:,2) = SM_dMRI_P2;
SM_P(:,3) = SM_dMRI_P3;
SM_P(:,4) = SM_dMRI_P4;
SM_P(:,5) = SM_dMRI_P5;

STD_P(:,1) = STD_dMRI_P1;
STD_P(:,2) = STD_dMRI_P2;
STD_P(:,3) = STD_dMRI_P3;
STD_P(:,4) = STD_dMRI_P4;
STD_P(:,5) = STD_dMRI_P5;

% Define experimental parameters
TE     = 51; % ms
Delta  = 35;
delta  = 9;

gyroMagnRatio =  267.5153151*10^(-6);

g = [166.8, 182.7, 197.3, 210.95, 235.85];
% g = [166.766843, 182.683923, 197.321189, 210.95, 235.85];
% b_vec =  gyroMagnRatio^2 * g.^2 * delta^2 * (Delta - delta/3); %ms/um^2
b_vec = [5, 6, 7, 8, 10];
q_vec = sqrt(b_vec./(Delta - delta/3))./delta;

for i=1:5 % For each ROI number
    display(['Analysis for ROI #' num2str(i)])
    Signal = SM_P(:,i);
    R_i    = R{i};
    
    % Estimate the effective radius from the measured (spherical mean)
    % dMRI signal using a pure diffusion model
    x0 = [1,    1];
    lb = [0.1,  0];
    ub = [10,  10000];
    [x,resnorm,residual,exitflag,output] = lsqcurvefit(@VanGelderenS_mod_vect_pure_diff, x0, b_vec(:), Signal(:), lb, ub);
    disp('------------------------------------------')
    display(['The effective radius estimated from the dMRI signal is ', num2str(x(1))])
    % Estimated from the exponential model
    Rmean_exp(i) = x(1);
    Signal_predicted(:,i) = VanGelderenS_mod_vect_pure_diff(x, b_vec(:));
    disp('------------------------------------------')
    
    if i==1
        Rxx = [R11; R12; R13; R14; R15];
    elseif i==2
        Rxx = [R11; R12; R13; R14; R15];
    elseif i==3
        Rxx = [R31; R32; R33; R34; R35];
    elseif i==4
        Rxx = [R41; R42; R43; R44; R45];
    elseif i==5
        Rxx = [R51; R52; R53; R54; R55];
    end
    
    % Generate the theoretical dMRI signal from a distribution of cylinder radii (for each ROI)
    S_theory = integrate_prop_signal(Rxx, b_vec(:));

    x0 = 1;
    lb = 0;
    ub = 100000;
    [scale,resnorm,residual,exitflag,output] = lsqcurvefit(@scale_factor, x0, S_theory(:), Signal(:), lb, ub);
    S_theory = S_theory * scale;
    S_theory_matrix(:,i) = S_theory(:);
    
    % Estimate the effective radius from the generated theoretical (spherical mean)
    % dMRI signal using a pure diffusion model
    x0 = [1,    1];
    lb = [0.1,  0];
    ub = [10,  10000];
    [x,resnorm,residual,exitflag,output] = lsqcurvefit(@VanGelderenS_mod_vect_pure_diff, x0, b_vec(:), S_theory(:), lb, ub);
    disp('------------------------------------------')
    display(['The (theoretical) effective radius is ', num2str(x(1))])
    Hmean_exp(i) = x(1);
    disp('------------------------------------------')
end

color(1,:) = [0 0.4470 0.7410];
color(2,:) = [0.9290 0.6940 0.1250];
color(3,:) = [0.4940 0.1840 0.5560];
color(4,:) = [0.4660 0.6740 0.1880];
color(5,:) = [0.6350 0.0780 0.1840];

%  ------------------------------------------------------------------------
figure('Renderer', 'painters', 'Position', [200 200 600 600])
% plot(Eff_32, Rmean_exp, 'o', 'color', 'b')
hold on
mdl = fitlm(Hmean_exp, Rmean_exp);                    % Fit Data
B   = mdl.Coefficients.Estimate;                      % Coefficients
CI  = coefCI(mdl);
% Coefficient Confidence Intervals
x = 0:0.1:6;
[Ypred,YCI] = predict(mdl, x');                      % Fitted Regression Line & Confidence Intervals
figure(1)

for i=1:length(Rmean_exp)
    h(i) = scatter(Hmean_exp(i), Rmean_exp(i), 200, 'p', 'MarkerFaceColor' , color(i,:), 'MarkerEdgeColor' , color(i,:));
end

hold on
pt = plot(x, Ypred,'-', x, YCI, '--r', 'LineWidth', 2, 'Color', [1, 0, 0, 0.5]);
hold off
grid
xlim([0, 6]); ylim([0, 6]);

hold on

axis('square')
plot(0:6,0:6, '-.k', 'LineWidth', 1)
xticks(1:6)
yticks(1:6)
xlabel('\it r_{eff-SEM-D} (\mu\itm)', 'FontSize',14)
ylabel('\it r_{eff-MRI-D} (\mu\itm)', 'FontSize',14)
title([' \it r_{eff-MRI-D}', ' ', 'vs' ' ', '\it r_{eff-SEM-D}'], 'fontsize', 16)

xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.92  * (yl(2)-yl(1)) + yl(1);
caption = sprintf('r_{eff-MRI-D} = %.2f + %.2f * r_{eff-SEM-D}', B(1), B(2));
text(xt, yt, caption, 'FontSize', 14, 'Color', 'red', 'FontAngle', 'italic');

xl = xlim;
yl = ylim;
xt = 0.06  * (xl(2)-xl(1)) + xl(1);
yt = 0.02  * (yl(2)-yl(1)) + yl(1);
str = texlabel('line: y=x');   
txt = text(xt, yt, str,  'FontSize', 14, 'FontAngle', 'italic');
set(txt,'Rotation', 45);

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.4         );

set( gca                       , ...
    'FontName'   , 'Helvetica' );

lgd = legend(h, {'Phantom 1','Phantom 2','Phantom 3', 'Phantom 4', 'Phantom 5', '', '', ''});
set(lgd,'Interpreter','latex', 'FontName', 'AvantGarde', 'Location', 'southeast');
lgd.FontSize = 14;

print('Figures/Figure6_A_pure_diffusion_mono_exp_rho_eq0', '-dpng', '-r600')

%  ------------------------------------------------------------------------

figure('Renderer', 'painters', 'Position', [200 200 600 600])
hold on
for i =1:5
    pt(i) = errorbar(b_vec*1000, log(SM_P(:,i)), STD_P(:,i), '.', 'color', color(i,:), 'MarkerSize', 20);
    plot(b_vec*1000, log(S_theory_matrix(:,i)),'--', 'color', color(i,:), 'LineWidth', 1.4);
end

axis('square')
xlabel('\it b-value (s/mm^2)', 'FontSize',14)
ylabel('\it log(S) (a.u)', 'FontSize',14)
title('dMRI signal vs SEM-based signal', 'fontsize', 16)

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

set(gca, ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.4         );

set( gca                       , ...
    'FontName'   , 'Helvetica' );

lgd = legend(pt, {'Phantom 1', 'Phantom 2', 'Phantom 3', 'Phantom 4', 'Phantom 5'});
set(lgd,'Interpreter','latex', 'FontName', 'AvantGarde', 'Location', 'northeast');
lgd.FontSize = 14;

set(lgd,'color','none');

title(lgd,'dMRI signal')

print('Figures/Figure6_B_logSignal_vs_b_value_mono_exp_rho0', '-dpng', '-r600')

%  ----------------      private function --------------------------------%
function S = integrate_prop_signal(Rxx, b_vec)

% I am summing over all signals generated for each radius in the sample
r2m   = sum(Rxx.^2);

S     = 0;
for i=1:length(Rxx)
    r = Rxx(i);
    signal_r2   =  VanGelderenS_mod_vect_pure_diff([r, 1], b_vec);
    % Pure diffusion model
    S  = S + (r^2/r2m) .* signal_r2; % standard theory
end
end

function Svg = scale_factor(x, S)
        scale = x;
        Svg = scale * S;
end
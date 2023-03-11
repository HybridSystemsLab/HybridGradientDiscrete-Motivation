%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clear all
close all
clc

%% Parameters

% plant parameters
theta = [1; 1];

% estimator parameters
s = 2*pi/30;               % sample period during flows
gammac = 0.5;
gammad = 100;

parameters.theta = theta;
parameters.s = s;
parameters.gammac = gammac;
parameters.gammad = gammad;


%% Create hybrid system
sys = Motivation(parameters);

tspan = [0, 10*pi];
jspan = [0, 500];
config = HybridSolverConfig('RelTol', 1e-4, 'AbsTol', 1e-4,  ...
                            'MaxStep', s/4, 'Refine', 4);

z0 = [0; 0];            % plant state
thetahat0 = [0;0];      % parameter estimate
tau_s0 = s;             % timer for samples during flows
tau0 = 0;               % timer for jumps
q0 = 0;                 % logic variable for jumps
thetahatc0 = thetahat0; % parameter estimate using only flows
thetahatd0 = thetahat0; % parameter estimate using only jumps
x0 = [z0;thetahat0;tau_s0;tau0;q0;thetahatc0;thetahatd0];
sol = sys.solve(x0, tspan, jspan, config);


%% Plots
HybridPlotBuilder.defaults.set('flow line width', 3, ...
                               'jump line width', 2,...
                               'label size', 34,...
                               'tick label size', 20,...
                               't_label', '$t$ [s]')
legendSize = 20;

figure; clf;
tlt = tiledlayout(2, 1);

nexttile(1)
hold on; grid on;
HybridPlotBuilder()....
    .flowColor('blue')...
    .jumpMarker('none')...
    .labels('$\psi_1$')...
    .tLabel('')...
    .plotFlows(sol.select(1))
set(gca,'Box','on');

nexttile(2)
hold on; grid on;
HybridPlotBuilder()....
    .flowColor('blue')...
    .jumpMarker('none')...
    .labels('$\psi_2$')...
    .plotFlows(sol.select(2))

xlim([0, 10*pi]) 
xticks([0 2*pi 4*pi 6*pi 8*pi 10*pi])
xticklabels({'$0$','$2\pi$','$4\pi$','$6\pi$','$8\pi$','$10\pi$'})
set(gca,'Box','on');

tlt.Padding = "none";
tlt.TileSpacing = "compact";
pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), 1.3*pos(4)])
movegui(gcf,'north')

%%
inD1 = (sol.x(:,6) >= 2*pi) | (sol.x(:,7) == 1);    % get points in the jump set
isFlowing = ~inD1;                                  % get points during flows
isJumping = inD1 + cat(1,0,inD1(1:end-1));          % get points either side of each jump

figure; clf;
hold on; grid on;

% plot thetahatc
HybridPlotBuilder().... 
    .flowColor('green')...% plot only the flows
    .jumpColor('none')...
    .plotFlows(sol.select(8:9),@(x) norm(x - theta))
HybridPlotBuilder().... % plot the jumps due to samples during flows
    .flowColor('none')...
    .jumpColor('green')...
    .jumpLineStyle('-')...
    .jumpMarker('none')...
    .filter(isFlowing)...
    .plotFlows(sol.select(8:9),@(x) norm(x - theta))
HybridPlotBuilder().... % plot the jumps due to the jump set
    .flowColor('none')...
    .jumpColor('red')...
    .jumpMarker('none')...
    .filter(isJumping)...
    .plotFlows(sol.select(8:9),@(x) norm(x - theta))

% plot thetahatd
HybridPlotBuilder().... % plot only the flows
    .flowColor('[0.9290 0.6940 0.1250]')...
    .jumpColor('none')...
    .plotFlows(sol.select(10:11),@(x) norm(x - theta))
HybridPlotBuilder().... % plot the jumps due to samples during flows
    .flowColor('none')...
    .jumpColor('[0.9290 0.6940 0.1250]')...
    .jumpLineStyle('-')...
    .jumpMarker('none')...
    .filter(isFlowing)...
    .plotFlows(sol.select(10:11),@(x) norm(x - theta))
HybridPlotBuilder().... % plot the jumps due to the jump set
    .flowColor('none')...
    .jumpColor('red')...
    .jumpMarker('none')...
    .filter(isJumping)...
    .plotFlows(sol.select(10:11),@(x) norm(x - theta))

% plot thetahat
HybridPlotBuilder().... % plot only the flows
    .flowColor('blue')...
    .jumpColor('none')...
    .labels('$|\tilde\theta_s|$')...
    .plotFlows(sol.select(3:4),@(x) norm(x - theta))
HybridPlotBuilder().... % plot the jumps due to samples during flows
    .flowColor('none')...
    .jumpColor('blue')...
    .jumpLineStyle('-')...
    .jumpMarker('none')...
    .labels('$|\tilde\theta_s|$')...
    .filter(isFlowing)...
    .plotFlows(sol.select(3:4),@(x) norm(x - theta))
HybridPlotBuilder().... % plot the jumps due to the jump set
    .flowColor('none')...
    .jumpColor('red')...
    .jumpMarker('none')...
    .labels('$|\tilde\theta_s|$')...
    .filter(isJumping)...
    .plotFlows(sol.select(3:4),@(x) norm(x - theta))

hpb = HybridPlotBuilder();
graphic_handles(1) = plot(nan, nan,'color',[0 1 0],'LineWidth',hpb.settings.flow_line_width);
graphic_handles(2) = plot(nan, nan,'color',[0.9290 0.6940 0.1250],'LineWidth',hpb.settings.flow_line_width);
graphic_handles(3) = plot(nan, nan,'color',[0 0 1],'LineWidth',hpb.settings.flow_line_width);
legend(graphic_handles,'Discretized continuous GD','Discrete GD','Discretized hybrid GD','Location','NorthEast','interpreter','latex','FontSize',legendSize)

ylim([0 1.7])
xlim([0, 10*pi]) 
xticks([0 2*pi 4*pi 6*pi 8*pi 10*pi])
xticklabels({'$0$','$2\pi$','$4\pi$','$6\pi$','$8\pi$','$10\pi$'})
set(gca,'Box','on');

tlt.Padding = "none";
tlt.TileSpacing = "compact";
pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), pos(4)])
movegui(gcf,'north')
set(gca, 'LooseInset', get(gca,'TightInset'))

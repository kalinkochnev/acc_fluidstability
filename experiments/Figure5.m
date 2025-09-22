% This file loads the maximum growth rate data and computes the
% eigendecomposition of P(t)
close all; clear; clc;
addpath("../utilities/")

txt = fileread('../Data/lyapunov_max_growth_n=2000.json');   % read entire file as text
loaded_debug = jsondecode(txt); 

omega = 0.5;

Pt_mats_cell = loaded_debug.P_opt;
n_timesteps = size(Pt_mats_cell, 1);
tvals = linspace(0, 2*pi/omega, n_timesteps);
n_vars = 5;

Pt_mats = zeros(n_vars, n_vars, n_timesteps);
for i = 1:n_timesteps
    Pt_mats(:, :, i) = squeeze(Pt_mats_cell(i, :, :));
end

% Eigen decomposition
[V, D] = eigenshuffle(Pt_mats);   % V: n_vars x n_vars x n_timesteps, D: n_vars x n_timesteps

% Extract dominant eigenvector at each time step
max_eigvectors = zeros(n_vars, n_timesteps);
for i = 1:n_timesteps
    [~, max_i] = max(D(:, i));
    max_eigvectors(:, i) = V(:, max_i, i);
end
% Get the moment in time where the eigenvalues are largest
[maximums, indices] = maxk(D(1, :), 2);
tpeak_1 = tvals(indices(1));
tpeak_2 = tvals(indices(2));
fprintf("The first peak for lambda1 occurs at t=" + tpeak_1 + " with magnitude " + maximums(1) + "\n\n");
fprintf("The second peak for lambda1 occurs at t=" + tpeak_2 + " with magnitude " + maximums(2) + "\n\n");
fprintf("The cosine minimum t=" + 2*pi/omega / 2 + '\n\n');
fprintf("The cosine maximum t=" + 2*pi/omega +'\n\n');


%------------------------------
% Plotting
%------------------------------
fig1 = figure;

lineStyles = {'b-', 'r--', 'k-', 'm--'};   % 'solid' → '-', '..' → ':'

% Subplot 1: Eigenvector components
set(gca, 'FontSize', 14);
grid on;
for r = 1:n_vars-1
    hold on;
    plot(linspace(0, 2*pi/omega, n_timesteps), max_eigvectors(r,:), lineStyles{r}, 'LineWidth', 1.2);
end
w_plot = plot(linspace(0, 2*pi/omega, n_timesteps), max_eigvectors(5,:), ':', 'LineWidth', 1.2);
w_plot.Color = '#00841a';
%w_plot.LineStyle = ':';
%w_plot.Marker = '';
lgd = legend({'$T$','$S$','$u$','$v$','$w$'}, 'Interpreter','latex', 'Location', 'northeastoutside');
pos = lgd.Position;     % [x y width height]
pos(1) = pos(1) - 0.25; % shift left by 5% of figure width
lgd.Position = pos;

xlabel('Time ($t$)','Interpreter','latex')
xticks(0:1:2*pi/omega);
xlim([0 2*pi/omega]);
ylabel('Eigenvector Component','Interpreter','latex')
print(fig1, 'Figure5.png', '-dpng', '-r300');

% Subplot 2: Eigenvalues
% Subplot 2: Eigenvalues
fig2 = figure;
set(gca, 'FontSize', 14);

% First eigenvalue line in black
plot(tvals, D(1,:), 'k-', 'LineWidth', 1.5);
hold on;
plot(tvals, D(2,:), '--', 'LineWidth', 1.5);

% Remaining eigenvalue lines
%for r = 2:n_vars
%    plot(tvals, D(r,:), lineStyles{r}, 'LineWidth', 1.2);
%end

xlabel('Time (t)','Interpreter','latex')
xticks(0:1:2*pi/omega);
xlim([0 2*pi/omega]);
ylabel('Eigenvalue','Interpreter','latex')
grid on;

% Add the background velocity with the second y-axis
fparams = experiment1_common();

yyaxis right
ax = gca;
set(ax, 'FontSize', 14);
ax.YColor = "blue";
plot(tvals, fparams.Au(tvals), '-.', 'LineWidth', 1.5, 'Color', ax.YColor);
ylabel('$A_U (t)$','Interpreter','latex')

%legend({'$\lambda_1$','','$\lambda_3$','$\lambda_4$','$\lambda_5$','$A_U(t)$'}, ...       'Interpreter','latex','Location','best')
legend({'$\lambda_1$', '$\lambda_2$', '$A_U(t)$'}, ...
       'Interpreter','latex','Location','best')


%print(fig2, 'Figure5.png', '-dpng', '-r300', '-loose');
clear all;
results.cvx_results = readmatrix("../Data/Convergence/lyapunov_convergence.csv");

txt = fileread('../Data/Convergence/params.json');   % read entire file as text
params = jsondecode(txt);
results.ode45_lambda = params.ode45_lambda;        % parse JSON into MATLAB struct/array
results.floquet_results = readmatrix("../Data/Convergence/floquet_convergence.csv");

% --- Full view ---
fig1 = figure;
ax = gca;
set(gca, 'FontSize', 25);

plot(results.cvx_results(:,1), results.cvx_results(:,2), ':', 'DisplayName', 'Convex Optimization', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(results.floquet_results(:,1), results.floquet_results(:,2), '--', 'DisplayName', 'Floquet', 'LineWidth', 1.5, 'MarkerSize', 6);
yline(results.ode45_lambda, '-', 'DisplayName', 'Numerical Simulation', 'Color', 'black', 'LineWidth', 2.5);
grid on;
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Growth Rate $$\lambda$$', 'Interpreter', 'latex', 'FontSize', 14);
ylim([0 1.2*results.ode45_lambda]);
legend('Location', 'southeast');
%set(fig1, 'outerposition', [1, 1, 1600, 1200]);
print(fig1, 'Figure4.png', '-dpng', '-r300');

% --- Zoomed view ---
fig2 = figure;
ax = gca;
set(gca, 'FontSize', 25);

plot(results.cvx_results(:,1), results.cvx_results(:,2), ':', 'DisplayName', 'Convex Optimization', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(results.floquet_results(:,1), results.floquet_results(:,2), '--', 'DisplayName', 'Floquet', 'LineWidth', 1.5, 'MarkerSize', 6);
yline(results.ode45_lambda, '-', 'DisplayName', 'Numerical Simulation', 'Color', 'red', 'LineWidth', 2.5);
grid on;
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Growth Rate $$\lambda$$', 'Interpreter', 'latex', 'FontSize', 14);
ylim([0 1.2*results.ode45_lambda]);
legend('Location', 'southeast');

xlim([min(results.cvx_results(:,1)), min(results.cvx_results(:,1)) + 100]); 
ylim([0, 1.2*results.ode45_lambda]);

%set(fig2, 'outerposition', [1, 1, 1600, 1200]);
%print(fig2, 'paper figures/convergence_zoom.png', '-dpng', '-r300');

paths = {
	"../Data/lyapunov_n=400.csv"
	"../Data/lyapunov_n=600.csv"
	"../Data/ode45 50 random initial conditions.csv"
	"../Data/floquet_n=2000.csv"
};
nFiles = numel(paths);

% X, Y coordinate limits
xlims = [-0.5 0.5];
ylims = [0 1.5];

x = linspace(xlims(1), xlims(2), 126);
y = linspace(ylims(1), ylims(2), 126);

% Preload
allData = cell(nFiles,1);
vmin = inf; vmax = -inf;

for i = 1:nFiles
    data = readmatrix(paths{i}); 
    data(data < -4) = NaN;
    allData{i} = data;


    vmin = min(vmin, min(data(:)));
    vmax = max(vmax, max(data(:)));
end
vmin=-4;
% Apply threshold to ode45 figure

% Create figure with tiled layout
f = figure('Units','inches','Position',[1 1 22 4.5]); % wide figure
t = tiledlayout(1, nFiles, 'TileSpacing','loose', 'Padding','loose');

for i = 1:nFiles
    nexttile;
    imagesc(x, y, allData{i});
    axis xy;
    %axis equal tight xy;
    clim([vmin vmax]);

    set(gca, 'FontSize', 14);
    xlabel("$k$", 'Interpreter','latex','FontSize', 24);
    ylabel("$m_0$", 'Interpreter','latex','FontSize', 24);
    xticks(linspace(xlims(1), xlims(2), 5));
    yticks(ylims(1):0.2:ylims(2));
	ylim([0 1.2]);
end

colormap("jet");
cb = colorbar;
cb.Layout.Tile = 'east';   % spans all plots
cb.Label.String = '$\log_{10}(\lambda)$'; % optional colorbar label
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 20;

% Compute the largest growth rates for all of these
N = 5; % number of maxima to extract

for i = 1:nFiles
    disp("Filename: " + paths{i});

    A = allData{i}(:);                % flatten matrix
    [maxVals, linearIdx] = maxk(A,N); % top N values + indices

    [rows, cols] = ind2sub(size(allData{i}), linearIdx);

    for j = 1:length(maxVals)
        fprintf("Max #%d at (k=%g, m0=%g) with lambda=%g\n", ...
                j, x(cols(j)), y(rows(j)), maxVals(j));
    end
    fprintf("\n");
end




% Save
savefig(f, "./paper figures/heatmap_all.fig");
print(f, "./paper figures/heatmap_all.png", "-dpng", "-r300");

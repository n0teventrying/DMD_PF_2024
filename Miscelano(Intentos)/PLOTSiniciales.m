D = 0.0001; % m^2/s
L_x = 2; % m
L_y = 1; % m
dx = 0.02; % m
dy = 0.02; % m

N_x = round(L_x / dx); % Number of nodes in x
N_y = round(L_y / dy); % Number of nodes in y

u_0 = zeros(N_y, N_x); % Note: N_y and N_x swapped to match Y and X meshgrid dimensions

[X, Y] = meshgrid((0:N_x-1) * dx, (0:N_y-1) * dy);
u_0 = exp(-1000 * ((X - L_x / 8).^2 + (Y - L_y / 2).^2));

figure;

% Plot using X and Y
pcolor(X, Y, u_0); 
shading interp;
colormap(flipud(hot));

% Set the axis limits to match the dimensions L_x and L_y
xlim([0 L_x]);
ylim([0 L_y]);

% Add filled contour lines
hold on;
contourf(X, Y, u_0, 10, 'LineColor', 'none');
hold off;

% Remove upper and right axis bars
box off;

% Label x and y axis
xlabel('X [m]');
ylabel('Y [m]');

% Add color bar and adjust its position
cbar = colorbar; % Assign colorbar handle for labeling
cbar.Label.String = 'Concentration'; % Labeling the colorbar
cbar.Label.Position = [cbar.Label.Position(1), max(Y(:))+0.03, cbar.Label.Position(3)]; % Adjust position above the color bar

% Add title
title('t = 0.000 [s]');







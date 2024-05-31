% PARAMETERS
D = 0.001; % m^2/s, diffusion coefficient
N_x = 100; % Number of nodes in the x-direction
N_y = 100; % Number of nodes in the y-direction
L_x = 2; % m, length of the domain in the x-direction
L_y = 3; % m, length of the domain in the y-direction
dx = L_x / N_x; % m, spatial step in the x-direction
dy = L_y / N_y; % m, spatial step in the y-direction
time = 20; % seconds, total simulation time
u_x = 0.00; % m/s, velocity in the x-direction
v_y = 0.05; % m/s, velocity in the y-direction
dt = (dx^2 * dy^2) / (2 * D * (dx^2 + dy^2)) / 10; % s, time step
t_nodes = floor(time / dt); % Number of time nodes

% Additional velocity field W parameters
w_x = 0.00; % m/s, additional velocity in the x-direction
w_y = 0.05; % m/s, additional velocity in the y-direction

% Initialize the contaminant concentration
u_0 = zeros(N_x, N_y);
[X, Y] = meshgrid((0:N_x-1) * dx, (0:N_y-1) * dy);
u = u_0;

% Simulation parameters
a = D * (dt / dx^2);
b = D * (dt / dy^2);
p = u_x * (dt / dx);
q = v_y * (dt / dy);

% Visualizing with a plot
figure;
h1 = pcolor((0:N_x-1) * dx, (0:N_y-1) * dy, u');
shading interp;
colormap(flipud(hot)); % Invert the colormap
colorbar_handle = colorbar;
caxis([0 1]);
hold on;
contour_lines = contour(X, Y, u', 10); % Initial contour plot

% Add a label to the colorbar
colorbar_handle.Label.String = 'C[mol/m^3]';

% New grid for W vector field
[W_X, W_Y] = meshgrid(linspace(0, L_x, 10), linspace(0, L_y, 10));
quiver(W_X, W_Y, w_x * ones(size(W_X)), w_y * ones(size(W_Y)), 'b'); % Quiver plot for W

xlim([0 L_x]);
ylim([0 L_y]);
xlabel('x [m]');
ylabel('y [m]');
title('Distribution at t: 0.000 [s].');
pause(0.01);

% Simulating
counter = 0;

while counter < time
    u_next = zeros(N_x, N_y);

    % Compute the next time step
    for i = 2:N_x-1
        for j = 2:N_y-1
            u_next(i, j) = u(i, j) + a * (u(i+1, j) - 2*u(i, j) + u(i-1, j)) + ...
                                      b * (u(i, j+1) - 2*u(i, j) + u(i, j-1)) - ...
                                      p * (u(i+1, j) - u(i, j)) - q * (u(i, j+1) - u(i, j));
        end
    end

    % Applying Dirichlet boundary conditions (fixed concentration)
    u_next(1, :) = 0;           % Left boundary
    u_next(end, :) = 0;         % Right boundary
    u_next(:, 1) = 1;           % Bottom boundary
    u_next(:, end) = 0;         % Top boundary

    % Update the concentration field
    u = u_next;

    counter = counter + dt;

    disp(['t: ', num2str(counter, '%.3f'), ' [s], Average concentration: ', num2str(mean(u(:)), '%.10f'), ' mol/m^3']);

    % Updating the plot
    set(h1, 'CData', u');
    caxis([0 max(u(:))]); % Adjust color limits to the maximum value of u

    % Update the contour plot
    if ishandle(contour_lines)
        delete(contour_lines); % Remove the old contour lines
    end
    contour_lines = contour(X, Y, u', 10); % Add new contour lines

    title(['Distribution at t: ', num2str(counter, '%.3f'), ' [s].']);
    %pause(0.01);
end

% Parameters
D = 0.001; % m^2/s
L_x = 2; % m
L_y = 1; % m
dx = 0.02; % m
dy = 0.02; % m
time = 15; % seconds
u_x = 0.07; % m/s
v_y = 0.00; % m/s
% Additional velocity field W parameters
w_x = 0.07; % m/s
w_y = 0.00; % m/s


N_x = round(L_x / dx); % Number of nodes in x
N_y = round(L_y / dy); % Number of nodes in y

%dt = min(dx^2 / (4 * D), dy^2 / (4 * D)); % s
%t_nodes = floor(time / dt);
% Time step and stability condition
dt = (dx^2 * dy^2) / (2 * D * (dx^2 + dy^2));
t_nodes = floor(time / dt);
% Initialize the contaminant concentration
[X, Y] = meshgrid((0:N_x-1) * dx, (0:N_y-1) * dy);

u_0 = exp(-1000 * ((X - L_x / 1.1).^2 + (Y - L_y / 2).^2))';
u = u_0;

% Visualization
figure('Position', [100, 100, 560, 420]); % Set the figure size to match the expected frame size
h = pcolor((0:N_x-1)*dx, (0:N_y-1)*dy, u'); % Use physical dimensions
shading interp;
colormap(flipud(hot));
colorbar;
caxis([0 1]);

hold on;
contour_lines = contour(X, Y, u', 10); % Initial contour plot
% Add a label to the colorbar
a = colorbar;
a.Label.String = 'C[mol/m^3]';

% New grid for W vector field
[W_X, W_Y] = meshgrid(linspace(0, L_x, 10), linspace(0, L_y, 10));
quiver(W_X, W_Y, w_x*ones(size(W_X)), w_y*ones(size(W_Y)), 'b'); % Quiver plot for W

xlim([0 L_x]);
ylim([0 L_y]);
xlabel('x[m]');
ylabel('y [m]')

title('Distribution at t: 0.000 [s].');



% Simulation parameters
a = D * (dt / dx^2);
b = D * (dt / dy^2);
p = u_x * (dt / dx);
q = v_y * (dt / dy);

% Simulation loop
counter = 0;

for t = 1:t_nodes
    u_next = zeros(size(u));
    
    % Periodic boundary conditions
    for i = 1:N_x
        for j = 1:N_y
            ip = mod(i, N_x) + 1;
            im = mod(i - 2, N_x) + 1;
            jp = mod(j, N_y) + 1;
            jm = mod(j - 2, N_y) + 1;
            
            u_next(i, j) = u(i, j) * (1 - 2*a - 2*b) ...
                + u(ip, j) * (a - p/2) + u(im, j) * (a + p/2) ...
                + u(i, jp) * (b - q/2) + u(i, jm) * (b + q/2);
        end
    end
    
    % Inyectando sustancia
    u = u_next + u_0 / 4;

    counter = counter + dt;

    disp(['t: ', num2str(counter, '%.3f'), ' [s], Average concentration: ', num2str(mean(u(:)), '%.10f'), ' mol/kg^3']);

    % Update plot
    set(h, 'CData', u');
    caxis([0 max(u(:))]); % Adjust color limits to the maximum value of u

    % Update the contour plot
    if ishandle(contour_lines)
        delete(contour_lines); % Remove the old contour lines
    end
    contour_lines = contour(X, Y, u', 10); % Add new contour lines
    title(['Distribution at t: ', num2str(counter, '%.3f'), ' [s].']);
    pause(0.1);
end


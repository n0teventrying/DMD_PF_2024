% Paramteros
D = 0.0001; % m^2/s
L_x = 2; % m
L_y = 1; % m
dx = 0.02; % m
dy = 0.02; % m
time = 80; % seconds
u_x = 0.007; % m/s
v_y = 0.00; % m/s

N_x = round(L_x / dx); % Number of nodes in x
N_y = round(L_y / dy); % Number of nodes in y

dt = min(dx^2 / (4 * D), dy^2 / (4 * D)); % s
t_nodes = floor(time / dt);

% Initialize the contaminant concentration
[X, Y] = meshgrid((0:N_x-1) * dx, (0:N_y-1) * dy);

  u_0 = exp(-1000 * ((X - L_x / 1.1).^2 + (Y - L_y /2).^2))';
u = u_0;

% Visualization
figure;
h = pcolor(u');
shading interp;
colormap(flipud(hot));
colorbar;
caxis([0 1]);
title('Distribution at t: 0.000 [s].');
pause(0.01);

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
            %TIPO 9 
            u_next(i, j) = u(i, j) * (1 - 2*a - 2*b) ...
                + u(ip, j) * (a - p/2) + u(im, j) * (a + p/2) ...
                + u(i, jp) * (b - q/2) + u(i, jm) * (b + q/2);
        end
    end
    
    % Inyectando sustancia
    u = u_next + u_0/6 ;

    counter = counter + dt;

    disp(['t: ', num2str(counter, '%.3f'), ' [s], Average concentration: ', num2str(mean(u(:)), '%.10f'), ' mol/kg^3']);

    % Update plot
    set(h, 'CData', u');
    title(['Distribution at t: ', num2str(counter, '%.3f'), ' [s].']);
    pause(0.1);
end

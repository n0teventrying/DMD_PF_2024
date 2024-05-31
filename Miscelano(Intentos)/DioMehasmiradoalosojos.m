%PARAMETROS
D = 0.001; % m^2/s

%Corrección, dar el número de nodos
N_x = 100;
N_y = 100;
L_x = 2; % m
L_y = 2; % m
dx = L_x/N_x; % m
dy = L_y /N_y; % m
time = 80; % seconds
u_x = 0.07; % m/s
v_y = 0.05; % m/s
%NÚMERO DE PUNTOS EN EL ARREGLO ESPACIAL
%N_x = round(L_x / dx); % Number of nodes in x
%N_y = round(L_y / dy); % Number of nodes in y
%NUMERO DE PUNTOS EN EL ARREGLO TEMPORAL Y CONDICIÓN DE ESTABILIDAD
%dt = 0.0004;
%dt = min(dx^2 / (4 * D), dy^2 / (4 * D))*0.001; % s
dt = (dx^2*dy^2)/(2*D*(dx^2+dy^2))/10;
t_nodes = floor(time / dt);
% CONCENTRACIÓN INICIAL DEL CONTAMINANTE
u_0 = zeros(N_x, N_y);
[X, Y] = meshgrid((0:N_x-1) * dx, (0:N_y-1) * dy);
u_0 = exp(-1000* ((X - L_x / 8).^2 + (Y - L_y / 2).^2))';
u = u_0;
% GRÁFICA
figure;
h = pcolor(u');
shading interp;
colormap (flipud(hot));
colorbar;
caxis([0 1]);
title('t: 0.000 [s].');
pause(0.01);
% DIFERENCIAS FINITAS
%constantes
a = D * (dt / dx^2);
b = D * (dt / dy^2);
p = u_x * (dt / dx);
q = v_y * (dt / dy);
% LOOP
counter = 0;
for t = 1:t_nodes
    u_next = zeros(size(u));
    % Tipo 9
    %Central en x; Central en y 
    u_next(2:end-1, 2:end-1) = u(2:end-1, 2:end-1) * (1 - 2*a - 2*b) + ...
        u(3:end, 2:end-1) * (a - p/2) + u(1:end-2, 2:end-1) * (a + p/2) + ...
        u(2:end-1, 3:end) * (b - q/2) + u(2:end-1, 1:end-2) * (b + q/2);
    %Tipo 1
    % Forward en x; Foward en y
    u_next(1, 1) = u(2, 1) * (-p - 2*a) + u(1, 2) * (-q - 2*b) + ...
        u(3, 1) * a + u(1, 3) * b + u(1, 1) * (1 + a + b + p + q);
    %Tipo 2
    % Backward in x and forward in y
    u_next(N_x, 1) = u(N_x, 1) * (1 - p + q + a + b) + ...
        u(N_x-1, 1) * (p - 2*a) + u(N_x, 2) * (-q - 2*b) + ...
        u(N_x-2, 1) * a + u(N_x, 3) * b;
    %Tipo 3
    % Backward en x; Backward en y
    u_next(N_x, N_y) = u(N_x, N_y) * (1 - p - q + a + b) + ...
        u(N_x-1, N_y) * (p - a) + u(N_x, N_y-1) * (q - 2*b) + ...
        u(N_x-2, N_y) * a + u(N_x, N_y-2) * b;
    %Tipo 4
    % Forward en x; Backward en y
    u_next(1, N_y) = u(1, N_y) * (1 + p - q + a + b) + ...
        u(2, N_y) * (-p - 2*a) + u(1, N_y-1) * (q - b) + ...
        u(3, N_y) * a + u(1, N_y-2) * b;
    %Tipo 5
    % Central en x; Forward en y
    u_next(2:end-1, 1) = u(2:end-1, 1) * (1 + q - a + b) + ...
        u(3:end, 1) * (a - p/2) + u(1:end-2, 1) * (a + p/2) + ...
        u(2:end-1, 2) * (-q - 2*b) + u(2:end-1, 3) * b;
    %Tipo 6
    % Backward en x; Central en y
    u_next(N_x, 2:end-1) = u(N_x, 2:end-1) * (1 - p + a - 2*b) + ...
        u(N_x-1, 2:end-1) * (p - 2*a) + u(N_x, 3:end) * (b - q/2) + ...
        u(N_x, 1:end-2) * (b + q/2) + u(N_x-2, 2:end-1) * a;
    %Tipo 7 
    % Central en x; Backward en y
    u_next(2:end-1, N_y) = u(2:end-1, N_y) * (1 - q - 2*a + b) + ...
        u(3:end, N_y) * (a - p/2) + u(1:end-2, N_y) * (a + p/2) + ...
        u(2:end-1, N_y-1) * (q - 2*b) + u(2:end-1, N_y-2) * b;
    %Tipo 8 
    % Forward en x, Central en y
    u_next(1, 2:end-1) = u(1, 2:end-1) * (1 + p + a - 2*b) + ...
        u(2, 2:end-1) * (-p - 2*a) + u(1, 3:end) * (b - q/2) + ...
        u(1, 1:end-2) * (b + q/2) + u(3, 2:end-1) * a;
    
    %Inyectando sustancia continuamente
    u = u_next +(u_0/4);

    counter = counter + dt;

    disp(['t: ', num2str(counter, '%.3f'), ' [s], Average concentration: ', num2str(mean(u(:)), '%.10f'), ' mol/kg^3']);

    % Updating the plot
    set(h, 'CData', u');
    title(['t: ', num2str(counter, '%.3f'), ' [s].']);
    pause(0.1);
end

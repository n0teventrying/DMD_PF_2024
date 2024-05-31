% Initialization
D = 0.0000; % m^2/s
L_x = 2; % m
L_y = 1; % m
dx = 0.02; % m
dy = 0.02; % m
time = 80; % seconds
u_x = 0.007; % m/s
v_y = 0.00; % m/s

N_x = round(L_x / dx); % Number of nodes in x
N_y = round(L_y / dy); % Number of nodes in y

dt = 0.0004;
%dt = min(dx^2 / (4*D ), dy^2 / (4 * D)); % s
t_nodes = floor(time / dt);

% Initialize the contaminant concentration
u_0 = zeros(N_x, N_y);
[X, Y] = meshgrid((0:N_x-1) * dx, (0:N_y-1) * dy);
u_0 = exp(-1000 * ((X - L_x / 8).^2 + (Y - L_y / 2).^2))';
u = u_0;

% Video setup
videoFileName = 'Advection.mp4';
v = VideoWriter(videoFileName, 'MPEG-4');
v.FrameRate = 10; % Adjust the frame rate as needed
open(v);

% Set up the figure
figure('Position', [100, 100, 560, 420]); % Set the figure size to match the expected frame size
h = pcolor((0:N_x-1)*dx, (0:N_y-1)*dy, u'); % Use physical dimensions
shading interp;
colormap(flipud(hot));
colorbar;
caxis([0 1]);
hold on;



% Simulation parameters
a = D * (dt / dx^2);
b = D * (dt / dy^2);
p = u_x * (dt / dx);
q = v_y * (dt / dy);

% Simulating
counter = 0;
capture_interval = 100; % Capture every 10th frame

for t = 1:t_nodes
    u_next = zeros(size(u));
    
    % Central nodes
    u_next(2:end-1, 2:end-1) = u(2:end-1, 2:end-1) * (1 - 2*a - 2*b) + ...
        u(3:end, 2:end-1) * (a - p/2) + u(1:end-2, 2:end-1) * (a + p/2) + ...
        u(2:end-1, 3:end) * (b - q/2) + u(2:end-1, 1:end-2) * (b + q/2);

    % Forward in x and y
    u_next(1, 1) = u(2, 1) * (-p - 2*a) + u(1, 2) * (-q + 2*b) + ...
        u(3, 1) * a + u(1, 3) * b + u(1, 1) * (1 + a + b + p + q);
    
    % Backward in x and forward in y
    u_next(N_x, 1) = u(N_x, 1) * (1 - p + q + a + b) + ...
        u(N_x-1, 1) * (p - 2*a) + u(N_x, 2) * (-q - 2*b) + ...
        u(N_x-2, 1) * a + u(N_x, 3) * b;
    
    % Backward in x and y
    u_next(N_x, N_y) = u(N_x, N_y) * (1 - p - q + a + b) + ...
        u(N_x-1, N_y) * (p - a) + u(N_x, N_y-1) * (q - 2*b) + ...
        u(N_x-2, N_y) * a + u(N_x, N_y-2) * b;
    
    % Forward in x and backward in y
    u_next(1, N_y) = u(1, N_y) * (1 + p - q + a + b) + ...
        u(2, N_y) * (-p - 2*a) + u(1, N_y-1) * (q - b) + ...
        u(3, N_y) * a + u(1, N_y-2) * b;
    
    % Central in x, forward in y
    u_next(2:end-1, 1) = u(2:end-1, 1) * (1 + q - a + b) + ...
        u(3:end, 1) * (a - p/2) + u(1:end-2, 1) * (a + p/2) + ...
        u(2:end-1, 2) * (-q - 2*b) + u(2:end-1, 3) * b;
    
    % Backward in x, central in y
    u_next(N_x, 2:end-1) = u(N_x, 2:end-1) * (1 - p + a - 2*b) + ...
        u(N_x-1, 2:end-1) * (p - 2*a) + u(N_x, 3:end) * (b - q/2) + ...
        u(N_x, 1:end-2) * (b + q/2) + u(N_x-2, 2:end-1) * a;
    
    % Central in x, backward in y
    u_next(2:end-1, N_y) = u(2:end-1, N_y) * (1 - q - 2*a + b) + ...
        u(3:end, N_y) * (a - p/2) + u(1:end-2, N_y) * (a + p/2) + ...
        u(2:end-1, N_y-1) * (q - 2*b) + u(2:end-1, N_y-2) * b;
    
    % Forward in x, central in y
    u_next(1, 2:end-1) = u(1, 2:end-1) * (1 + p + a - 2*b) + ...
        u(2, 2:end-1) * (-p - 2*a) + u(1, 3:end) * (b - q/2) + ...
        u(1, 1:end-2) * (b + q/2) + u(3, 2:end-1) * a;
    
    % Update u by adding u_next and the constant source u_0
    u = u_next ;

    counter = counter + dt;

    disp(['t: ', num2str(counter, '%.3f'), ' [s], Average concentration: ', num2str(mean(u(:)), '%.10f'), ' mol/kg^3']);

    % Updating the plot
    set(h, 'CData', u');
    title([' t: ', num2str(counter, '%.3f'), ' [s].']);
    


    % Capture every nth frame
    if mod(t, capture_interval) == 0
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    
    pause(0.01);
end
    
% Capture the plot as a frame in the video
%    frame = getframe(gcf);
 %   writeVideo(v, frame);
    
  %  pause(0.1);
%end

% Close the video file
close(v);

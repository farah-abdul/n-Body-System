%% NAME : ANGELA BENZIGAR
%% PURPOSE : 3-BODY SIMULATION 

% PROVIDED PARAMETERS
vx = 0.347112814; 
vy = 0.532726852; 
T = 6.325; % Period of orbit
dt = T / (1 * 10^5); % Calculated time step
timeVector = 0:dt:T; % Calculated time vector

% PROVIDED CONDITIONS
x0 = [-1 0 1]; 
y0 = [0 0 0]; 
u0 = [vx -2*vx vx]; 
v0 = [vy -2*vy vy]; 

% RUN THE SIMULATION USING EXPLICIT EULER FUNCTION 
[x,y] = n_body_sim(timeVector, x0, y0, u0, v0); 

% PLOT THE RESULTS
figure(1)
hold on 
plot(x(:,1), y(:,1), 'r-', 'LineWidth', 1);
hold on
plot(x(:,2), y(:,2), 'g--', 'LineWidth', 1);
hold on
plot(x(:,3), y(:,3), 'b:', 'LineWidth', 1); 
hold on
xlabel('x position, x(t)') 
ylabel('y position, y(t)'); 
title("Three Body Problem Simulation using Explicit Euler Method")
legend('Body 1', 'Body 2', 'Body 3'); 
hold off
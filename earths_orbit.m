% need it to output time t, pos x, pos y, xvelocity u, yvelocity v, x acceleration, y acceleration
% solving for x acc and y acc uses x,y, so acceleration doesnt need to be an output
% vy0 (initial y velocity) and Nt (# of timesteps) will be inputs

function [t, x, y, u, v] = earths_orbit(v0, Nt)

%% initial conditions and parameters

P = 31558150;                % Total Period, sec (t0 = 0, tf = 351558150 --> not sure if I should add that to my script)
dt = P/Nt;                   % delta t (time step), sec
Msun = 1.989e30;             % Mass of Sun, kg
G = 6.674e-11;               % Gravitational Constant, m^3/kg/s^2
C = G*Msun;                  % Constant
x0 = 1.471e11;               % initial x position, m
y0 = 0;                      % initial y position, m
u0 = 0;                      % initial x velocity, m/s
% like i said, v0 is an input

%% initialize matrices so that MATLAB doesnt get angry at me later

x = zeros(Nt+1, 1);          % will get Nt+1 values of x
y = zeros(Nt+1, 1);          % will get Nt+1 values for y
u = zeros(Nt+1, 1);          % will get Nt+1 values of u
v = zeros(Nt+1, 1);          % will get Nt+1 values for v
t = [0:Nt]'*dt;              % time at n in each iteration, transposed to make column vector
x(1) = x0;                 
y(1)= y0;                    
u(1) = u0;
v(1) = v0;

for n = 1:Nt

    r = sqrt( x(n)^2 + y(n)^2 ); % radius

    au = -C * x(n)/r^3;         % x acceleration, m/s^2
    av = -C * y(n)/r^3;         % y acceleation, m/s^2

    % half step x, y, u, ,v (xh = xn + dt/2*K1)
    xh = x(n) + (dt/2)*u(n);
    yh = y(n) + (dt/2)*v(n);
    uh = u(n) + (dt/2)*au;
    vh = v(n) + (dt/2)*av;

    % half step k values, first two K values found from the last step (because u and v are technically K values for x and y)

    rh = sqrt( xh^2 + yh^2 );           % half step radius, m
    auh = -C*(xh/rh^3);                 % half step x acceleration, m/s^2
    avh = -C*(yh/rh^3);                 % half step y acceleration, m/s^2

    % full step x, y, u, v (xh = xn + dt/2*K2)

    x(n+1) = x(n) + dt*uh;
    y(n+1) = y(n) + dt*vh;
    u(n+1) = u(n) + dt*auh;
    v(n+1) = v(n) + dt*avh;
end
end


   
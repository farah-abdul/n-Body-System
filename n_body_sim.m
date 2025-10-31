function [x_hist, y_hist] = n_body_sim(t, x0, y0, u0, v0, substeps)

%check how many bodies and default to no substep if < 6
if nargin < 6
    substeps = 1;
end
% define variables
n = numel(x0);
Nt = numel(t);
dt = t(2) - t(1);
% initialize vectors
x_hist = zeros(Nt,n);
y_hist = zeros(Nt,n);

x = x0(:).';  % row vector
y = y0(:).';  % row vector
u = u0(:).';  % row vector
v = v0(:).';  % row vector

x_hist(1,:) = x;
y_hist(1,:) = y;

% Acceleration vector calculation function
    function [ax,ay] = accel(x,y)
        dx = x.' - x;
        dy = y.' - y;
        r2 = dx.^2 + dy.^2;
        r2(1:n+1:end) = Inf; % don't apply to own body
        invr3 = 1 ./ (r2 .* sqrt(r2));
        ax = -sum(dx .* invr3, 2).';
        ay = -sum(dy .* invr3, 2).';
    end

[ax, ay] = accel(x, y);

% Explicit Euler with substepping
for k = 2:Nt
    dt_sub = dt / substeps;
    %Run substeps
    for s = 1:substeps
        % half step velocities
        u = u + 0.5*dt_sub*ax;
        v = v + 0.5*dt_sub*ay;

        % full step positions
        x = x + dt_sub*u;
        y = y + dt_sub*v;

        % update accelerations
        [ax, ay] = accel(x, y);

        % half step velocities again
        u = u + 0.5*dt_sub*ax;
        v = v + 0.5*dt_sub*ay;
    end

    % store new positions
    x_hist(k,:) = x;
    y_hist(k,:) = y;
end
end

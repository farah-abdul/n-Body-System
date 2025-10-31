%% Aidan Reckamp part g)
clear
clc
n_values = [2, 4, 8, 16, 32];
vcirc = zeros(size(n_values));

Nt = 20000;
T_sim = 0.2*2*pi;
T_sim_initial = 0.01*2*pi;

for k = 1:length(n_values)
    n = n_values(k);

    %initial positions
    theta = linspace(0, 2*pi, n+1); theta(end) = [];
    x0 = cos(theta);
    y0 = sin(theta);

    %initial guess
    %Special Case for n=32 please work
    if n == 32
        v0_guess = 4.316239;
        max_iter = 150;
        substeps = 30;
        Nt_sim = 20000;
        T_seg_initial = 0.002*2*pi;
        T_seg_long = 0.02*2*pi;
        vel_factor = 0.6;
    else
        %regular case
        v0_guess = sqrt(n);
        max_iter = 100;
        substeps = 5;
        Nt_sim = Nt;
        T_seg_initial = T_sim_initial;
        T_seg_long = T_sim;
        vel_factor = 0.3;
    end
    %go to max iteration and create new function inputs
    for iter = 1:max_iter
        u0 = -v0_guess * y0;
        v0_vec = v0_guess * x0;

        % make first iteration short
        if iter == 1
            T_seg = T_seg_initial;
        else
            T_seg = T_seg_long;
        end
        t = linspace(0, T_seg, Nt_sim);

        % simulate wtih function
        [x_hist, y_hist] = n_body_sim(t, x0, y0, u0, v0_vec, substeps);

        % radial deviation from unit circle calculation
        r_all = sqrt(x_hist.^2 + y_hist.^2);

        % update velocity guess
        scale = 1/sqrt(mean(r_all(:).^2));
        v0_guess = v0_guess * (1 + vel_factor*(scale-1)); %scale factor for n=32 case

        % print catch
        fprintf('n=%d, Iter %d: v_guess=%.6f, max radial dev=%.6f\n', ...
            n, iter, v0_guess, max(abs(r_all(:)-1)));

        % end if radial deviation is within tolerance
        if max(abs(r_all(:)-1)) < 1e-5
            break;
        end
    end
    % print final value
    vcirc(k) = v0_guess;
    fprintf('n = %d, v_circ = %.5f, max radial dev = %.5f\n', ...
        n, vcirc(k), max(abs(r_all(:)-1)));

    % plot trajectories to check
    figure;
    hold on;
    for j = 1:n
        plot(x_hist(:,j), y_hist(:,j), '-');
        plot(x_hist(1,j), y_hist(1,j), 'o', 'MarkerFaceColor','k');
    end
    axis equal
    xlabel('x')
    ylabel('y')
    title(sprintf('Particle trajectories for n = %d (v_{circ} = %.4f)', n, vcirc(k)))
    grid on;
    hold off;
end

% plot vcirc vs n figure 
figure;
plot(n_values, vcirc, 'o-', 'LineWidth', 1.5);
xlabel('Number of particles n');
ylabel('Circular orbit velocity v_{circ}');
title('Circular-orbit velocity vs number of particles');
grid on;

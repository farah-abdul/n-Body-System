% need to find optimal v0, so we're using earths_orbit

%% Given/Known

G = 6.674e-11;                            % Gravitational constant, m^3/kg/s^2
Msun = 1.989e30;                          % Mass of Sun, kg
C = G*Msun;                               % constant

% velocity will be highest at perihelion and lowest at aphelion due to gravitational pull
% we calculate vperi and vaph to get our sample range for v, and give a +- 2km margin of error

rperi = 147.1e9;                          % distance to sun at perihelion (closest point during orbit), m
raph = 152.1e9;                           % distance to sun at aphelion (furthest point during orbit), m

% using ellipse speed, not circular, for more accurate result, because we care about physics here

a = (C*(P/(2*pi))^2)^(1/3);               % semimajor axis, m
vperi = sqrt(C*((2/rperi)-(1/a)));        % ellipse speed at perihelion, m/s
v_upper = vperi + 2000;                   % vperi with a +2km margin (upper bound), m/s 
vaph = sqrt(C*((2/raph)-(1/a)));          % ellipse speed at aphelion, m/s
v_lower = vaph - 2000;                    % vaph with a -2km margin (lower bound), m/s
vrange = v_lower:v_upper;                 % range of values we're checking
tol = 10^-5;                              % need at least 0.0001 m/s precision              
Nt_vector = [ 10 100 1000 10000 100000 ];
P = 31558150;                             % Period, s
dt = P./ Nt_vector;                       % delta t, s
opt = zeros(length(Nt_vector),3);         % initializing an array that saves Nt, optimal velocity, and minimum distance after each outer loop iteration
                                          
%% OUTER LOOP: looping each Nt value

for k =1:length(Nt_vector)                % repeat for every Nt option
    Nt = Nt_vector(k);                    % pulling one Nt value from vector of Nt values

    % implementing a coarse search in order to narrow down our search radius
    % coarse search

    v_coarse = v_lower:10:v_upper;        % coarse search to get a smaller range to work with

    % initialize coarse search

    v_opt = v_coarse(1);                                    % initializing a "best candidate" before looping (starting from first iteration), m/s
    [t,x,y,u,v] = earths_orbit(v_opt, Nt);                 
    
    % setting this as the minimum distance and comparing the rest of the iterations against it
    % if they're better, then they replace it, and the corresponding v_opt replaces v0 to be used in the next iteration
    min_dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);     % (final x- initial x), (final y - initial y), gives distance between initial and final positions using the current optimal velocity

    % now that we have a baseline for min distance and optimal v, we create a loop to calculate the same parameters, but for every velocity within our range and step size
    % will compare results to the previous iteration

    % create loop that cycles through the v_coarse options and compares each iteration to the last

    %% INNER LOOP: looping through every v0 candidate to find the optimal one for given Nt
    for n = 2:length(v_coarse)
        v0 = v_coarse(n);
        [t,x,y,u,v] = earths_orbit(v0, Nt);
        dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
        % time to check if this is better than the previous iteration
        if dist < min_dist                                  % want the distance between initial and final position to be as small as possible
            min_dist = dist;                                % updates min_dist estimates if current iteration shows smaller error
            v_opt = v0;                                     % updates v0 to the velocity corresponding to the smaller error (smaller min_dist)
        end
    end

    % more precise loop to achieve a result within our tolerance

    v_precise = (v_opt-50):(v_opt+50);                      % small range around our best candidate, small increments
    for n = 1:length(v_precise)
        v0 = v_precise(n);
        [t,x,y,u,v] = earths_orbit(v0, Nt);
        dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
        if dist < min_dist
            min_dist = dist;                                % updates min_dist estimates if current iteration shows smaller error
            v_opt = v0;                                     % updates v0 to the velocity corresponding to the smaller error (smaller min_dist)
        end
    end

    % create an even MORE precise loop to get a result within our tolerance 10e-5
    % step size selected to get precise results without extending runtime

    v_tiny = (v_opt-1):0.1:(v_opt+1);                       % small range around our best fine candidate, small increments
    for n = 1:length(v_tiny)
        v0 = v_tiny(n);
        [t,x,y,u,v] = earths_orbit(v0, Nt);
        dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
        if dist < min_dist
            min_dist = dist;                                % updates min_dist estimates if current iteration shows smaller error
            v_opt = v0;                                     % updates v0 to the velocity corresponding to the smaller error (smaller min_dist)
        end  
    end

    % almost at our tolerance

    v_super_tiny = (v_opt-0.1):0.001:(v_opt+0.1);
    for n = 1:length(v_super_tiny)
        v0 = v_super_tiny(n);
        [t,x,y,u,v] = earths_orbit(v0, Nt);
        dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
        if dist < min_dist
            min_dist = dist;                                % updates min_dist estimates if current iteration shows smaller error
            v_opt = v0;                                     % updates v0 to the velocity corresponding to the smaller error (smaller min_dist)
        end
    end

    % finally at our tolerance 10^-5 ! 

     v_extremely_tiny = (v_opt-0.001):0.0001:(v_opt+0.001);
    for n = 1:length(v_extremely_tiny)
        v0 = v_extremely_tiny(n);
        [t,x,y,u,v] = earths_orbit(v0, Nt);
        dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
        if dist < min_dist
            min_dist = dist;                                % updates min_dist estimates if current iteration shows smaller error
            v_opt = v0;                                     % updates v0 to the velocity corresponding to the smaller error (smaller min_dist)
        end
    end

    % plotting at this point shows a decrease in position error up until Nt=10^4 (error ~140m ), then it jumps up at Nt=10^5 (error ~2900m )
    % kind of funky, no?
    % calculations show the difference between optimal velocity at Nt = 10^4 and Nt = 10^5 occurs at the 4th to 5th decimal place
    % MATLAB recognizes the results as exactly the same due to this. It sees a smaller dt lead to the same results as the next largest dt
    % So Nt = 10^5 appears to have a larger error than it really does 
    % This means search is still (somehow) too coarse
 

    % getting too ambitious, more precise loop to catch difference between Nt = 10^4 and Nt = 10^5

    v_ambitious = (v_opt-0.0001):0.00001:(v_opt+0.0001); 
    for n = 1:length(v_ambitious)
        v0 = v_ambitious(n);
        [t,x,y,u,v] = earths_orbit(v0, Nt);
        dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
        if dist < min_dist
            min_dist = dist;                                % updates min_dist estimates if current iteration shows smaller error
            v_opt = v0;                                     % updates v0 to the velocity corresponding to the smaller error (smaller min_dist)
        end
    end

    % this loop catches the minute difference between Nt = 10^4 and Nt = 10^5 and yields a small error of ~16m for Nt = 10^5
    % could (and should definitely stop here)
    % but the line isn't as straight as I would like

    % enough is enough please stop

    v_pleasestop = (v_opt-0.00001):0.000001:(v_opt+0.00001); 
    for n = 1:length(v_pleasestop)
        v0 = v_pleasestop(n);
        [t,x,y,u,v] = earths_orbit(v0, Nt);
        dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
        if dist < min_dist
            min_dist = dist;                                % updates min_dist estimates if current iteration shows smaller error
            v_opt = v0;                                     % updates v0 to the velocity corresponding to the smaller error (smaller min_dist)
        end
    end

    % this loop reduced the Nt=10^4 error by a factor of 10, but the Nt=10^5 is still ~16m
    % guess we need another loop to get it down to near 0
    % least efficient code ever


    % okay seriously this is the last one

     v_LastTimeISwear = (v_opt-0.000001):0.0000001:(v_opt+0.000001);
    for n = 1:length(v_LastTimeISwear)
        v0 = v_LastTimeISwear(n);
        [t,x,y,u,v] = earths_orbit(v0, Nt);
        dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
        if dist < min_dist
            min_dist = dist;                                % updates min_dist estimates if current iteration shows smaller error
            v_opt = v0;                                     % updates v0 to the velocity corresponding to the smaller error (smaller min_dist)
        end
    end                             

    % Nt=10^5 error down to ~2m !!
    % we're close to 0 ...might as well ?


    % i really think i can get the error close to 0 so we're doing one more loop

     v_TooCarriedAway = (v_opt-0.0000001):0.00000001:(v_opt+0.0000001);
    for n = 1:length(v_TooCarriedAway)
        v0 = v_TooCarriedAway(n);
        [t,x,y,u,v] = earths_orbit(v0, Nt);
        dist = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2);
        if dist < min_dist
            min_dist = dist;                                % updates min_dist estimates if current iteration shows smaller error
            v_opt = v0;                                     % updates v0 to the velocity corresponding to the smaller error (smaller min_dist)
        end
    end

    % graph shows an almost straight line of steadily decreasing error
    % Nt=10^5 error is practically 0
    % can relax now

    opt(k,:) = [Nt, v_opt, min_dist];                       % will save the Nt, v_opt, and min_dist found at current Nt's loop
end

 % this is just v_opt, but I created a new variable so that it doesnt overwrite the other iterations
 v_optimal = opt(end,2);                                    % opt(end,2) means last row (where Nt=100k) and second column (which holds v_opt values)

 % absolute difference between each v_opt (for each Nt) and the best v_opt across all Nt values
 abs_diff = abs(v_optimal - opt(:,2));                      
 min_dist_vector = opt(:,3);

%% PLOTS

figure(1)
% plotting the error between optimal velocity at each Nt value and overall optimal velocity across ALL Nt values (NOT the true, analytical velocity)
grid on;
% using log because Nt = 10, 100, 1000 etc
loglog(Nt_vector, abs_diff);                                % difference between v_opt and v_best for each Nt, against Nt
xlim([10, 10e4]);
ylim([1e-3, 1e4]);
xlabel('Nt');
ylabel('Deviation from Optimal Velocity (m/s)');    
title('Convergence to Optimal Velocity');

% difference should decrease because we're increasing Nt (aka decreasing time step) so it should be getting more accurate
% more samples = more accurate = more optimal results !

figure(2)
grid on;
loglog(Nt_vector, min_dist_vector);                         % minimum distance from each Nt, against the time step that produced each results
xlabel('Nt');
ylabel('Deviation from Initial Position (m)');
title('Position Error vs Time Step');

% How many initial samples did we place at a tolerance of 10^-5?
samples_tol = length(v_coarse) + length(v_precise) + length(v_tiny) + length(v_super_tiny) + length(v_extremely_tiny)

% How many did we place to get to the tolerance I wanted (because I got carried away)?
samples_total = length(v_coarse) + length(v_precise) + length(v_tiny) + length(v_super_tiny) + length(v_extremely_tiny) + length(v_ambitious) + length(v_pleasestop) + length(v_LastTimeISwear) + length(v_TooCarriedAway)

%% so what exactly is the posititon error showing us? Realistically, is Earth supposed to end up at the same position
% (since its 1 full year), and what we're doing is showing how we are only approximating its position at the end of the
% year with RK2? And how as we increase the # of time steps we can get a better approximation, which is why we're looking
% for the minimum distance between the initial and final position, because whatever is closest to 0 is closest to the
% true value, which is 0, because Earth will end up in the same position. I just don't understand why we had to cycle
% through different initial velocities for this
% we know the true value is at perihelion, so i guess the point of this is to show how with more time steps, our RK2
% converges towards the true value for the optimal initial velocity
% which just means the velocity the Earth travels to be able to return to its initial position in exactly 1 year
% error decreases as Nt goes up 
% MAGNITUDE OF ERROR? in best case scenario

% why binary search method not appropriate
% conceptually kind of similar (narrowing our range until we're satisfied)
% we are not necessarily looking for a root, and even if we were
% there could be false sign changes, leading us to discard the other stuff where the true soln migh t be
% too much variance with our funtion over the period with varying v values
% for v's that are too small, (physically) that would mean Earth would be closer to the Sun than initial position
% so the distance between final and initial increases
% for v too large, Earth escapes the gravitational field of the Sun and flies off into space
% once again, distance between initial and final would be too large
% minimum distance is somewhere in between
% with a binary method, we would see the minimum approaching closer to the center
% but we wouldn't know which side to discard. Will it still continue to go down for a few more terms?
% which would mean we havent found the minimum yet
% or is it going to start going up again?
% we wouldnt be able to know which piece to discard without doing extra computational work, defeating the entire purpose





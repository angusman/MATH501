function chainsim(N)
% CHAINSIM Simulate a free falling chain as a system of n-pendulums.
% Input: N is the number of links in chain.
%  If N = 1, single pendulum demo is run.
%  If N = 2, double pendulum demo is run.
%  For N > 2, solver stops when last bob is vertical.

%% Implementation
close all; shg
% Number of pendulums
global n g
n = N; % number of links
g = 1;  % gravity
 
% Initial Settings
A0 = makeA0(n);  % make the initial angle
v0 = zeros(n,1); % Angular velocity
y0 = [A0; v0];

if n == 1 || n == 2
    % These are demos to show that the solver works as expected
    tspan = [0 100];
    opts = odeset('maxstep',.02);
    [T,Y] = ode23(@odefun,tspan,y0,opts);
    step = 5; % info for plotter
else
    % We will detect when the last angle is vertical, otherwise the results
    % become unstable for large t...
    tspan = [0 inf];
    opts = odeset('maxstep',.02,'events',@eventsg);
    [T,Y] = ode23(@odefun,tspan,y0,opts);
    step = 1; % info for plotter
end

% computed angles
theta = Y(:,1:n);
% setting (x,y) postions cumultaively
x = cumsum(cos(theta-pi/2), 2); 
y = cumsum(sin(theta-pi/2), 2); 

% show the pendulum system
figure(1)
for t = 1:step:length(T)
    plot(x(1:t,n),y(1:t,n),'-') % trail of last bob
    hold on
    plot([0 x(t,1:n)],[0 y(t,1:n)],'r-o') % visual pendulum
    hold off
    axis(1.25*[-n n -n n])
    axis square
    title(['t = ' num2str(T(t))])
    pause(.01)
end

% Our second figure plots the height of the last bob vs. freefall
figure(2)
hold on
plot(T,y(1:length(T),n))
ylabel('height')
xlabel('time')

y0 = y(1,n);
plot(T,y0-0.5*g*T.^2,'r--')

legend('chain','freefall')




       


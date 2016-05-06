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
    opts = odeset('maxstep',.02,'events',@events);
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

function dy = odefun(t,y)
global n g
% We set up a system of the form M*dy = f.

% LHS (Constructing M)
I = eye(n);
M = zeros(2*n);

% Upper-left block
M(1:n,1:n) = I;

% Lower-right block
if n == 1
    M(2,2) = 1;
else
    for i = 1:n
        for k = 1:n
            if k <= i-1
                M(n+i,n+k) = cos(y(k)-y(i));
            elseif k == i
                M(n+i,n+k) = n-i+1;
            else % k >= i+1
                M(n+i,n+k) = (n-k+1)*cos(y(i)-y(k));
            end
        end
    end
end

% RHS (Constructing f)
f = zeros(2*n,1);

% First n entries
f(1:n) = y(n+1:2*n);

% Last n entries
if n == 1
    f(2) = -g*sin(y(1));
else
    for i = 1:n
        Sless = zeros(1,i-1); Smore = zeros(1,n-i);
        for j = 1:i-1
           Sless(j) = (0*2*y(n+j)*y(n+i) - y(n+j)^2)*sin(y(j)-y(i)); 
        end
        for j = i+1:n
           Smore(j-i) = (n-j+1)*sin(y(i)-y(j))*y(n+j)^2; 
        end
        f(n+i) = -g*(n-i+1)*sin(y(i)) - sum(Sless) - sum(Smore);
    end
end

% Solve
dy = M\f;

function A0 = makeA0(n)
    A0 = zeros(n,1);
    if rem(n,2) == 0
        fst = 1:n/2-1;
        mid = n/2;
        nnd = n/2+1:n;
    else
        fst = 1:(n-1)/2;
        mid = (n+1)/2;
        nnd = (n+3)/2:n;
    end
    
    % Rectangular Option
        A0(fst) = zeros(length(fst),1);
        A0(mid) = pi/2;
        A0(nnd) = pi*ones(length(nnd),1) - pi/20 * (length(nnd):-1:1)';
    
% Arc Option
% A0(1:end) = (pi/n)*(0:n-1)';

function [value,isterminal,direction] = events(t,y)
global n
value = y(n);
isterminal = 1;
direction = 0;
       


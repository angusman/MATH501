function [T,U] = npend
close all; shg
global n m L grav
%% n = 1
n = 3;
m = 1;
L = 1;
grav = 1;

% link numbers
% N = n-1;
% theta = pi*(0:1/N:1)';

theta = 0 * ones((n-1)/2,1);
theta = [theta; pi/2; pi/2+pi/2*(1/((n-1)/2):1/((n-1)/2):1)'];
thetadot = zeros(n,1);

u0 = [theta; thetadot];
  
opts = odeset('outputfcn',@swingplot,'maxstep',.02);
[T,U] = ode23(@odefun,[0 10],u0);

try
    ode23(@odefun,[0 10],u0,opts);
catch
    return
end


function du = odefun(t,u)
global n m L grav
% Mass Matrix
M = zeros(2*n);
M(1:n,1:n) = eye(n);

for i = 1:n
   for j = 1:n 
      M(i+n,j+n) = L^2*G(i,j)*cos(u(i)-u(j)); 
      if j == i
         M(i+n,i+n) = M(i+n,i+n) - (1/4)*L^2; 
      end
   end
end
% cond(M)

% RHS Vector
f = zeros(2*n,1);
f(1:n) = u(n+1:2*n);

for i = 1:n
    for j = 1:n
        f(n+i) = f(n+i) - L^2*G(i,j)*sin(u(j)-u(i))*u(j+n);
    end
    f(n+i) = f(n+i) - grav*sin(u(i))*G(i,i);
end

du = M\f;



function status = swingplot(t,u,task)
% Plot function for classic double pendulum.
global n grav
persistent plt plt2 titl orbt erasemode
% Coordinates of both bobs

theta = u(1:n);
x = cumsum(sin(theta));
y = cumsum(-cos(theta));

x2 = n;
y2 = -grav*t(end)^2;

hold on

switch task
    
    case 'init'
        
        % Initialize plot
        
        erasemode = verLessThan('matlab','8.4');
        
        % For x-y pal
        plt = plot([0; x],[0; y],'o-');
        plt2 = plot([x2],[y2],'ro');
        
        axis(n*2.25*[-1 1 -1 1])
        axis square
        if erasemode
            titl = title(sprintf('t = %8.1f',t),'erasemode','xor');
        else
            titl = title(sprintf('t = %8.1f',t));
        end
        
    case ''
        
        % Update plot
        
        set(plt,'xdata',[0; x],'ydata',[0; y])
        set(plt2,'xdata',[x2],'ydata',[y2])
        % Display time in title
        
        set(titl,'string',sprintf('t = %8.1f',t));
        drawnow
end

% Terminate ode solver after mouse click (not really)

status = ~isempty(get(gca,'userdata'));
pause(.01)
hold off


function g = G(i,k)
global n
g = n - (max(i,k) - 1/2);


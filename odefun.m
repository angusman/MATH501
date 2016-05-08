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


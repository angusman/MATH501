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
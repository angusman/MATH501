function [value,isterminal,direction] = eventsg(t,y)
global n
value = y(n);
isterminal = 1;
direction = 0;
end

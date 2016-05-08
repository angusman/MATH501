function acclerationgraph(nrange)
% finding acceleration for various numbers of chain links.
% as the chain sim is kinda inefficent expect this to take awhile
% Inputs:
%  nrange - range of numbers of chain links i.e. 5:15
x = [];
y = [];
for k = 1:range(nrange)+1
    x = [k+min(nrange)-1 x];
    y = [chainsim(k+min(nrange)-1,0) y];

end

figure(4)
hold on
title('Acceleration over chain length with g = 1')
plot(x,y,'o-r')
xlabel('links in chain')
ylabel('accleration of end bob')
hold off
    
end


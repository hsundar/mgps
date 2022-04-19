n = [68921, 185193, 531441,1442897,4173281,11390625];
setup = [0.68, 2.3, 4.1, 15.9, 33.1, 135];
solve1 = [0.13, 0.22, 1.08, 1.7, 8.7, 15.6];
solve2 = [0.06, 0.06, 0.3, 0.6, 2.2, 49];

plot(n,setup); hold on;
plot(n,solve1);
plot(n,solve2);

legend('setup', 'merge', 'leaf')

% for i=1:5
%   s = setup(i+1)/setup(i); m = solve1(i+1)/solve1(i); l = solve2(i+1)/solve2(i);
%   nn = n(i+1)/n(i);
%   fprintf('n: %g, s: %g, m: %g, l: %g\n', nn, s, m, l);
% end

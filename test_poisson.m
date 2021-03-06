k=6;
rhs = @(x,y,z) -3*k*k*pi^2*sin(k*pi*x).*sin(k*pi*y).*sin(k*pi*z);
test_op;
num_elems = [16,16,16];
order = 5;
mu = @(x,y,z)(1); % @water;

%-----
m = mgps.mesh(num_elems, @mgps.xform.identity);
t1 = tic;
m.initialize_leaves(order, op, rhs, mu);
t2 = toc(t1);
fprintf('Initialized leaves in %g secs\n', t2);

%num_bdy = prod(num_elems.*order + 1) - prod(num_elems)*(order-1)^3; 


num_elems  = prod(m.nelems);
nnf = m.refel.nnf;
num_bdy = nnf * ( 3*num_elems + m.nelems(1)*m.nelems(2) + m.nelems(2)*m.nelems(3) + m.nelems(1)*m.nelems(3));

d = m.trace_diagonal();
invdiag = 1./d;

% plot(d);
% 
% u = rand(num_bdy,1);
% r = m.trace_residual(u);
% %norm(r)
% 
% stem(r)

% jacobi test 

omega = 6/7;
u = zeros(num_bdy,1);
t1 = tic;
for i=1:20
  res  = invdiag .* m.trace_residual(u);
  u = u + omega.*res;
  r = norm(res);
  norm(r)
  %plot(u);
end
t2 = toc(t1);
fprintf('Solved %g secs\n', t2);

m.clear();
v = m.solve_leaf(u);
t1 = tic;
v = m.solve_leaf(u);
t2 = toc(t1);
fprintf('Solved leaf: %g secs\n', t2);

r = norm(res);
norm(r)

nlin = m.nelems(1)*order +1;
w = reshape(v, [nlin, nlin, nlin]);
imagesc(w(:,:,(nlin-1)/2)); colorbar;
%slice(w, (nlin-3)/2, (nlin-3)/2, (nlin-3)/2); colorbar; 

% function mu = water(x,y,z)
%   r2 = (x.*x + y.*y + z.*z);
%   if ( r2 < 0.36)
%     mu = 80;
%   elseif (r2 > 0.49)
%     mu = 4;
%   else
%     mu = 80 - (80.0 - 4.0)/(0.49 - 0.36)*(r2 - 0.36);
%   end
% end 

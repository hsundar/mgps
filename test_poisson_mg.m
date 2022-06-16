%% Problem specification  
% discretization
dim = 3; 
nelems = [8,16];   % [2,4,8,16]
orders = [3];      % 

% operator 
op = struct();
op.dxx = -1; op.dyy = -1; op.dzz = -1;
op.dxy = 0; op.dyz = 0; op.dxz = 0;
op.dx  = 0; op.dy  = 0; op.dz  = 0;
op.b = 0;

k = 3;
k2 = 7;
rhs = @(x,y,z) -3*k*k*pi*pi*sin(k*pi*x).*sin(k*pi*y).*sin(k*pi*z); % + -3*k2*k2*pi*pi*sin(k2*pi*x).*sin(k2*pi*y).*sin(k2*pi*z);
mu = @(x,y,z)(1);
bdy = @(x,y,z)(0);

xform = @mgps.xform.identity;
%---

%% Setup 

num_hgrids = length(nelems);
num_pgrids = length(orders);

num_grids = num_hgrids + num_pgrids - 1;

disp(['Creating h-grid: ' num2str(1) ' of ' num2str(num_grids) ', order = ' num2str(orders(1)) ', nelem = ' num2str(nelems(1))]);
m = mgps.mesh(repmat(nelems(1), 1, dim), xform);
coarse = mgps.grid(m, orders(1));

cgrid = coarse;

for i=2:num_hgrids
  disp(['Creating h-grid: ' num2str(i) ' of ' num2str(num_grids) ', order = ' num2str(orders(1)) ', nelem = ' num2str(nelems(i))]);
  m = mgps.mesh(repmat(nelems(i), 1, dim), xform);
  grid = mgps.grid(m, orders(1), coarse);
  coarse = grid;
end

hfine = nelems(num_hgrids);

% disp('Creating p-grids now');
for i=2:num_pgrids
  disp(['Creating p-grid: ' num2str(i+num_hgrids-1) ' of ' num2str(num_grids) ', order = ' num2str(orders(i)) ', nelem = ' num2str(hfine)]);
  m = mgps.mesh(repmat(hfine, 1, dim), xform);
  grid = mgps.grid(m, orders(i), coarse);
  coarse = grid;
end

% Now setup operators
grid.assemble_operators(op, mu, rhs, bdy);

grid.is_finest = true;
grid.Fine = [];



%% Solve FAS

% u = cgrid.get_u0();
% while ( ~isempty(cgrid) )
% %   disp('-------------------------------------');
% %   fprintf('  Level %d \n', cgrid.level);
% %   disp('-------------------------------------');
%   if ( ~isempty(cgrid.Coarse) )
%     u = cgrid.prolong( u );
%   end
%   res = cgrid.residual(cgrid.get_u0());
%   r1 = norm(res);
%   u1 = cgrid.smooth(300, cgrid.get_u0());
%   res = cgrid.residual(u1);
%   r2 = norm(res);
% 
%   res = cgrid.residual(u);
%   r3 = norm(res);
%   % [u, rr, iter] = cgrid.solve(20, 'jacobi', 3, 3, u);
%   u = cgrid.smooth(300, u);
%   res = cgrid.residual(u);
%   r4 = norm(res);
% 
%   u = u1;
%   fprintf(' Level %d norms, [%g, %g], [%g, %g]\n', cgrid.level, r1, r2, r3, r4);
%   cgrid = cgrid.Fine;
% end
% 
% res = grid.residual(u);
% norm(res)
% v = grid.solve_leaf(u);
% grid.plot_skel(u);

%% Solve 
% 

%grid.Coarse.pfac = .1;
grid.pfac = 20;
grid.set_smoother('jacobi');
% 
u = grid.get_u0();
%u = grid.smooth(120, u);
[u, rr, iter] = grid.solve(10, 'jacobi', 3, 3, grid.get_u0(), u);

% pfr = zeros(100,1);
% for pf=1451:1470
%   u = grid.get_u0();
%   grid.pfac = pf;
%   [u, rr, iter] = grid.solve(20, 'jacobi', 3, 3, u);
%   pfr(pf-1450) = rr;
% end
% figure,plot(1451:1470, pfr);
% [m,i] = min(pfr)
res = grid.residual(u, grid.get_u0());
norm(res)
% v = grid.solve_leaf(u);
% grid.plot_skel(u);
% max(v(:))
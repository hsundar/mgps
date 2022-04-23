%% Problem specification  
% discretization
dim = 3; 
nelems = [8,16];
orders = [5];

% operator 
op = struct();
op.dxx = -1; op.dyy = -1; op.dzz = -1;
op.dxy = 0; op.dyz = 0; op.dxz = 0;
op.dx  = 0; op.dy  = 0; op.dz  = 0;
op.b = 0;

k=2;
rhs = @(x,y,z) -12*pi^2*sin(k*pi*x).*sin(k*pi*y).*sin(k*pi*z);
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

%% Solve 

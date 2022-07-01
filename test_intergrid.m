%% Problem specification  
% discretization
dim = 3; 
nelems = [4, 8];
orders = [5];

% operator 
op = struct();
op.dxx = -1; op.dyy = -1; op.dzz = -1;
op.dxy = 0; op.dyz = 0; op.dxz = 0;
op.dx  = 0; op.dy  = 0; op.dz  = 0;
op.b = 0;

k=3;
rhs = @(x,y,z) -3*k*k*pi*pi*sin(k*pi*x).*sin(k*pi*y).*sin(k*pi*z);
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

%% Intergrid

u = grid.get_u0();

num_elem_c = grid.Coarse.Mesh.nelems;
num_elem_f = grid.Mesh.nelems;
nnf = grid.Mesh.refel.nnf;

forder = [1 3 2 5 4 6];

for k=1:num_elem_f(3)
  for j=1:num_elem_f(2)
    for i=1:num_elem_f(1)
      e = sub2ind (num_elem_f, i, j, k);
      fid = grid.Mesh.get_global_faces(e);
      for ff=1:6
        u(((fid(ff)-1)*nnf+1):(fid(ff)*nnf)) = forder(ff);
      end
    end
  end
end
grid.plot_skel(u);
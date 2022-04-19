classdef mesh < handle
  % MESH A container class for homg regular grid meshes
  
  properties (SetAccess = private)
    dim=2;
    nelems%=[8 8];
    order
    
    coords       % element vertices
    Xf           % transform
    
    % problem specific
    refel
    coeff
    muvec % maybe ?
    rhs % maybe ?
    % elementwise Solution and D2N operators.
    S;
    D2N;
  end % properties

  methods
    function mesh = mesh(nelems, X)
      mesh.dim    = length(nelems);
      mesh.nelems  = nelems;
      mesh.Xf = X;
      
      mesh.coeff = @(x,y,z)(1);
    end

    function plot(self)
      figure(1);
      c = [31/256,171/256,226/256];   % default color of grid
      lw = 1;                         % default line width of grid
      
      if (self.dim == 2 )
        [x,y] = ndgrid(0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0);
        pts = [x(:) y(:)];
        coords = self.Xf(pts);
        plot(coords(:,1), coords(:,2), 'ko');
        hold on;
        x = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        plot(x,y, 'Color',c,'LineWidth',lw);
        plot(x',y', 'Color',c,'LineWidth',lw);
        % boundary test
        idx = self.get_boundary_node_indices(1);
        plot(coords(idx,1), coords(idx,2), 'ro');
        title(['Quad Mesh - ' num2str(self.nelems(1)) 'x' num2str(self.nelems(2)) ])
        axis square
      else
        % 2 xy planes
        [x,y] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0 );
        
        % z = 0
        z = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(2)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        % z = 1
        z = ones  (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(2)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        %--------------------------------------------------------------------
        % 2 yz planes
        [y,z] = ndgrid( 0:1/self.nelems(2):1.0, 0:1/self.nelems(3):1.0 );
        % x = 0
        x = zeros (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(2)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(2)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(2)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        % x = 1
        x = ones (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(2)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(2)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(2)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        %--------------------------------------------------------------------
        % 2 xz planes
        [x,z] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(3):1.0 );
        % y = 0
        y = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        % y = 1
        y = ones (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)) ...
            );
        hold on;
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        
        % pretty views etc
        view(3); axis equal;
        title(['Hex Mesh ' num2str(self.nelems(1)) 'x' num2str(self.nelems(2)) 'x' num2str(self.nelems(3))])
      end
      set (gcf, 'renderer', 'opengl');
      cameratoolbar('show');
      cameratoolbar('setmode', 'orbit');
    end

    function plot_fx(self, fx)
      % display the mesh. Needs X.
      hFig = figure(1);
      set(gcf,'PaperPositionMode','auto')
      set(hFig, 'Position', [200 200 800 800])
      
      if (self.dim == 2 )
        scl=8;
        [x,y] = ndgrid(0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(2)):1.0);
        pts = [x(:) y(:)];
        coords = self.Xf(pts);
        
        ci = arrayfun( fx, coords(:,1), coords(:,2) );
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            zeros(size(x)), reshape(ci, size(x)), ...
            'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [x,y] = ndgrid(0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0);
        pts = [x(:) y(:)];
        coords = self.Xf(pts);
        x = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        plot(x,y, 'k-');
        plot(x',y', 'k-');
        
        axis square;
        view(0,90);
        colorbar;
      else
        % will draw 6 planes ...
        %--------------------------------------------------------------------
        % 2 xy planes
        scl=8;
        [x,y] = ndgrid( 0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(2)):1.0 );
        % z = 0
        z = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)), ...
            reshape(ci, size(x)), ...
            'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [x,y] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0 );
        z = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(2)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        % z = 1
        [x,y] = ndgrid( 0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(2)):1.0 );
        z = ones  (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)), ...
            reshape(ci, size(x)), ...
            'EdgeColor','none','LineStyle','none' ...
            ); % 'FaceColor', 'interp', 'FaceLighting', 'phong'
        hold on;
        [x,y] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(2):1.0 );
        z = ones  (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(2)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(2)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(2)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        %--------------------------------------------------------------------
        % 2 yz planes
        [y,z] = ndgrid( 0:1/(scl*self.nelems(2)):1.0, 0:1/(scl*self.nelems(3)):1.0 );
        % x = 0
        x = zeros (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)), ...
            reshape(ci, size(x)), ...
            'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [y,z] = ndgrid( 0:1/self.nelems(2):1.0, 0:1/self.nelems(3):1.0 );
        % x = 0
        x = zeros (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(2)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(2)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(2)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        % x = 1
        [y,z] = ndgrid( 0:1/(scl*self.nelems(2)):1.0, 0:1/(scl*self.nelems(3)):1.0 );
        x = ones (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)), ...
            reshape(ci, size(x)), ...
            'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [y,z] = ndgrid( 0:1/self.nelems(2):1.0, 0:1/self.nelems(3):1.0 );
        x = ones (size(y));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(2)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(2)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(2)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        %--------------------------------------------------------------------
        % 2 xz planes
        [x,z] = ndgrid( 0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(3)):1.0 );
        % y = 0
        y = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)), ...
            reshape(ci, size(x)), ...
            'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [x,z] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(3):1.0 );
        y = zeros (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        % y = 1
        [x,z] = ndgrid( 0:1/(scl*self.nelems(1)):1.0, 0:1/(scl*self.nelems(3)):1.0 );
        y = ones (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        ci = arrayfun( fx, coords(:,1), coords(:,2), coords(:,3) );
        surf(reshape(coords(:,1), size(x)), ...
            reshape(coords(:,2), size(x)), ...
            reshape(coords(:,3), size(x)), ...
            reshape(ci, size(x)), ...
            'EdgeColor','none','LineStyle','none' ...
            );
        hold on;
        [x,z] = ndgrid( 0:1/self.nelems(1):1.0, 0:1/self.nelems(3):1.0 );
        y = ones (size(x));
        pts = [x(:) y(:) z(:)];
        coords = self.Xf(pts);
        x1 = reshape(coords(:,1), self.nelems(1)+1, self.nelems(3)+1);
        y1 = reshape(coords(:,2), self.nelems(1)+1, self.nelems(3)+1);
        z1 = reshape(coords(:,3), self.nelems(1)+1, self.nelems(3)+1);
        plot3(x1,y1,z1, 'k-');
        plot3(x1',y1',z1', 'k-');
        
        % pretty views etc
        view(3); %axis square
        %view(150,40);
        colorbar;
      end
      
      % print ('-depsc2', 'warped-2d.eps');
      % matlab2tikz (['fan-' num2str(self.dim) 'd.tikz'], 'checkForUpdates', false, 'showInfo', false);
      set (gcf, 'renderer', 'opengl');
      cameratoolbar('show');
      cameratoolbar('setmode', 'orbit');
    end

    function plot_field(self, u)
      % function plot_field(self, u)
      %    will plot a solution on the grid

    end

    function plot_trace(self, lambda)
      % function plot_lambda(self, lambda)
      %    will plot a trace on the grid
      
    end

    function initialize_leaves(self, order, op, rhs)
      self.order = order;
      r = mgps.refel ( self.dim, order );
      self.refel = r;
            
      num_elems  = prod(self.nelems);
      
      % loop over elements
      for e=1:num_elems
        % idx = self.get_node_indices (e, order);  % needed ?
        pts = self.element_nodes(e, r);
        detJac = prod(0.5./self.nelems); %self.geometric_factors(refel, pts);

        %% setup RHS
        if ( isnumeric(rhs) && isscalar(rhs) )
          rhs = repmat(rhs, r.ni, 1);
        end

        if ( isa(rhs, 'function_handle') )
          rhs_val = feval(rhs, pts(r.ii,1), pts(r.ii,2), pts(r.ii,3));
        else
          rhs_val = rhs;
        end

        %% setup op
        A = zeros(r.n);

        for name = fieldnames(op).'
          if ( isa(op.(name{1}), 'function_handle') )
            A = A + feval(op.(name{1}), pts(:,1), pts(:,2), pts(:,3)) .* r.(name{1})(jac);
          elseif ( op.(name{1}) ~= 0 )
            A = A + op.(name{1})(:) .*  r.(name{1})(detJac);
          end
        end

        %% compute S and D2N maps
        dA = decomposition(A(r.ii, r.ii));
        Ainv = @(u) dA \ u;
        S = Ainv([-A(r.ii, r.ee), rhs_val(:)]);

        % Append boundary points to solution operator:
        tmpS = zeros(r.n, size(S, 2));
        tmpS(r.ii,:) = S;
        tmpS(r.ee,:) = eye(r.ne, r.ne+1);
        self.S{e} = tmpS(:,[r.f1s r.f3s r.f2s r.f5s r.f4s r.f6s end]);

        % D2N 
        D2N = zeros(r.ne, r.ne+1);
        D2N(r.f1s,:) = -r.Dx(r.f1, :) * self.S{e};
        D2N(r.f2s,:) = r.Dx(r.f2, :) * self.S{e};
        D2N(r.f3s,:) = -r.Dy(r.f3, :) * self.S{e};
        D2N(r.f4s,:) = r.Dy(r.f4, :) * self.S{e};
        D2N(r.f5s,:) = -r.Dz(r.f5, :) * self.S{e};
        D2N(r.f6s,:) = r.Dz(r.f6, :) * self.S{e};

        % nnf = length(r.f6s);
        % disp(nnf)
        self.D2N{e} = D2N; %mat2cell(D2N, [nnf nnf nnf nnf nnf nnf], [nnf nnf nnf nnf nnf nnf 1]);
      end % element loop
    end % initialize leaves


    function d = trace_diagonal(mesh)
      num_elems  = prod(mesh.nelems);
      nnf = mesh.refel.nnf;
      ndof = nnf * ( 3*num_elems + mesh.nelems(1)*mesh.nelems(2) + mesh.nelems(2)*mesh.nelems(3) + mesh.nelems(1)*mesh.nelems(3));
      d = zeros(ndof,1);
      for e=1:num_elems
        bdy_idx = [];
        fid = mesh.get_global_faces(e);
        for f=1:length(fid)
          bdy_idx = [bdy_idx, ((fid(f)-1)*nnf+1):(fid(f)*nnf)];
        end
        if (max(bdy_idx) > ndof)
          disp('exceeding bdy index');
        end
        d(bdy_idx) = d(bdy_idx) + diag(mesh.D2N{e});
      end
    end

    function r = trace_residual(mesh, u)
      % function r = trace_residual(mesh, u)
      num_elems  = prod(mesh.nelems);
      
      nnf = mesh.refel.nnf;

      r = zeros(size(u));
      % loop over elements
      for e=1:num_elems
        bdy_idx = [];
        % copy trace to element boundaries 
        fid = mesh.get_global_faces(e);
        for f=1:length(fid)
          bdy_idx = [bdy_idx, ((fid(f)-1)*nnf+1):(fid(f)*nnf)];
        end
        %    compute D2N*u 
        re = mesh.D2N{e}*[u(bdy_idx); 1];
        r(bdy_idx) = r(bdy_idx) - re;
      end % elem loop   
    end % trace_residual

    function u = solve_leaf(mesh, trace)
      num_elems  = prod(mesh.nelems);
      ndofs = prod(mesh.nelems.*mesh.order + 1);
      u = zeros(ndofs,1);
      nnf = mesh.refel.nnf;
      
      for e=1:num_elems
        idx = mesh.get_node_indices(e);
        bdy_idx = [];
        % copy trace to element boundaries 
        fid = mesh.get_global_faces(e);
        for f=1:length(fid)
          bdy_idx = [bdy_idx, ((fid(f)-1)*nnf+1):(fid(f)*nnf)];
        end
        u(idx) = mesh.S{e}*[trace(bdy_idx); 1];
      end
    end

    function clear(mesh)
      mesh.D2N = 0;
    end

    %%===================================================================

    function bdy_idx = get_element_boundary_node_indices(self, eid)
      % determine global node indices for a given element
      order = self.order;

      if ( self.dim == 2)
        [i,j] = ind2sub (self.nelems, eid);
       
        i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low   = (j-1)*order + 1;   j_high =  j*order + 1;
          
        [i,j] = ndgrid(i_low:i_high, j_low:j_high);
          
        idx     = sub2ind (self.nelems*order + 1, i(:), j(:));
      else
        [i,j,k] = ind2sub (self.nelems, eid);
           
        i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
        j_low   = (j-1)*order + 1;   j_high =  j*order + 1;
        k_low   = (k-1)*order + 1;   k_high =  k*order + 1;
            
        [i,j,k] = ndgrid(i_low:i_high, j_low:j_high, k_low:k_high);
            
        idx     = sub2ind (self.nelems*order + 1, i(:), j(:), k(:) );
      end

      bdy_idx = idx(self.refel.ee);
    end



    function [J, D] = geometric_factors( self, refel, pts )
      % Np =  refel.Nrp ^ mesh.dim;
      
      % compute x,y,z for element
      % pts = self.element_nodes(elem, refel);
      
      % change to using Qx etc ?
      if (refel.dim == 1)
          xr  = refel.Dg*pts;
          J = xr;
      elseif (refel.dim == 2)
          [xr, xs] = homg.tensor.grad2 (refel.Dg, pts(:,1));
          [yr, ys] = homg.tensor.grad2 (refel.Dg, pts(:,2));
          
          J = -xs.*yr + xr.*ys;
      else
          [xr, xs, xt] = homg.tensor.grad3 (refel.Dg, pts(:,1));
          [yr, ys, yt] = homg.tensor.grad3 (refel.Dg, pts(:,2));
          [zr, zs, zt] = homg.tensor.grad3 (refel.Dg, pts(:,3));
          
          J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
      end
      
      if (nargout > 1)
          if (refel.dim == 1)
              D.rx = 1./J;
          elseif (refel.dim == 2)
              D.rx =  ys./J;
              D.sx = -yr./J;
              D.ry = -xs./J;
              D.sy =  xr./J;
          else
              D.rx =  (ys.*zt - zs.*yt)./J;
              D.ry = -(xs.*zt - zs.*xt)./J;
              D.rz =  (xs.*yt - ys.*xt)./J;
              
              D.sx = -(yr.*zt - zr.*yt)./J;
              D.sy =  (xr.*zt - zr.*xt)./J;
              D.sz = -(xr.*yt - yr.*xt)./J;
              
              D.tx =  (yr.*zs - zr.*ys)./J;
              D.ty = -(xr.*zs - zr.*xs)./J;
              D.tz =  (xr.*ys - yr.*xs)./J;
          end
          
      end
  end

    % get boundary face number (no- shared edge/vertex nodes)
    function fid = get_global_faces(self, elem)
      % fid = get_global_faces(self, elem) 
      % 
      % ordering for element f1 f3 f2 f5 f4 f6 
      if (self.dim == 2)
        [i,j] = ind2sub (self.nelems, elem);
        % order all x=c edges followed by y=c edges
        f1 = (j-1)*(self.nelems(1)) + i;
        f2 = (j-1)*(self.nelems(1)) + i + 1;
        offset = self.nelems(2)*(self.nelems(1)+1);
        f3 = offset + i + (j-1)*(self.nelems(1));
        f4 = offset + i + j*(self.nelems(1));

        fid = [f1 f3 f2 f4];
      else
        [i,j,k] = ind2sub (self.nelems, elem);
        % all yz faces (x=c), followed by xz (y=c) and then xy (z=c) faces
        f1 = (k-1)*self.nelems(2)*(self.nelems(1)+1) + (j-1)*(self.nelems(1)+1) + i;
        f2 = (k-1)*self.nelems(2)*(self.nelems(1)+1) + (j-1)*(self.nelems(1)+1) + i + 1;
        %
        offset = self.nelems(3)*self.nelems(2)*(self.nelems(1)+1) ;
        f3 = offset + (k-1)*(self.nelems(2)+1)*self.nelems(1) + (j-1)*(self.nelems(1)) + i;
        f4 = offset + (k-1)*(self.nelems(2)+1)*self.nelems(1) + j*(self.nelems(1)) + i;
        %
        offset = offset +  self.nelems(3)*(self.nelems(2)+1)*self.nelems(1) ;
        f5 = offset + (k-1)*self.nelems(2)*self.nelems(1) + (j-1)*(self.nelems(1)) + i;
        f6 = offset + k*self.nelems(2)*self.nelems(1) + (j-1)*(self.nelems(1)) + i;
        %
        fid = [f1 f3 f2 f5 f4 f6];
      end

    end

    function idx = get_node_indices ( self, eid )
      % determine global node indices for a given element
      if ( self.dim == 2)
          [i,j] = ind2sub (self.nelems, eid);
          
          i_low   = (i-1)*self.order + 1;   i_high =  i*self.order + 1;
          j_low   = (j-1)*self.order + 1;   j_high =  j*self.order + 1;
          
          [i,j] = ndgrid(i_low:i_high, j_low:j_high);
          
          idx     = sub2ind (self.nelems*self.order + 1, i(:), j(:));
      else
          [i,j,k] = ind2sub (self.nelems, eid);
          
          i_low   = (i-1)*self.order + 1;   i_high =  i*self.order + 1;
          j_low   = (j-1)*self.order + 1;   j_high =  j*self.order + 1;
          k_low   = (k-1)*self.order + 1;   k_high =  k*self.order + 1;
          
          [i,j,k] = ndgrid(i_low:i_high, j_low:j_high, k_low:k_high);
          
          idx     = sub2ind (self.nelems*self.order + 1, i(:), j(:), k(:) );
      end
    end

    % todo: update for chebyshev
    function coords = element_nodes(self, elem, refel)
      h = 1./self.nelems;
      
      if ( self.dim == 2)
          [i,j] = ind2sub (self.nelems, elem);
          idx = [i j];
      else
          [i,j,k] = ind2sub (self.nelems, elem);
          idx = [i j k];
      end
      
      p_mid = (idx - 0.5) .* h;
      p_gll = refel.r * 0.5 * h;
      nodes = bsxfun(@plus, p_mid, p_gll) ;
      
      if ( self.dim == 2)
          [x, y] = ndgrid(nodes(:,1), nodes(:,2));
          pts = [x(:) y(:)];
      else
          [x, y, z] = ndgrid(nodes(:,1), nodes(:,2), nodes(:,3));
          pts = [x(:) y(:) z(:)];
      end
      
      coords = self.Xf(pts);
    end

    function idx = get_boundary_node_indices(self, order)
      % function idx = get_boundary_node_indices(self, order)
      %    returns indices of boundary nodes, for setting
      %    boundary conditions
      if (self.dim == 2)
        [x,y] = ndgrid(1:self.nelems(1)*order+1,1:self.nelems(2)*order+1);
          
        idx = [ find(x == 1);     find(x == (self.nelems(1)*order+1));
              find(y == 1);     find(y == (self.nelems(2)*order+1))  ];
          
        idx = unique(idx);
      else
        [x,y,z] = ndgrid(1:self.nelems(1)*order+1,1:self.nelems(2)*order+1,1:self.nelems(3)*order+1);
          
        idx = [ find(x == 1);     find(x == (self.nelems(1)*order+1));
              find(y == 1);     find(y == (self.nelems(2)*order+1));
              find(z == 1);     find(z == (self.nelems(3)*order+1))  ];
          
        idx = unique(idx);
      end
    end % get_boundary_node_indices 
  end % methods 
end % class
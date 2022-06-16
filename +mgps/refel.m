classdef refel < handle
  % REFEL hexahedral reference element
  % currently using Chebyshev
    
  properties
    dim
    N     % polynomial order
    p     % number of 1D interpolation points (N+1)
        
    r     % 1D reference coordinates of the chebyshev nodes ( 1 x p )
        
    D     % 1D Derivative 
          % ( p x p )
          % D(i,j) = chebyshev_i' (r_j)
    
    I
    D2
    Dx
    Dy
    Dz
    Dxx
    Dxy
    Dyy
    Dyz
    Dxz
    Dzz

    L

    %% indices
    n
    ni
    ne
    nnf
    ii % interior
    ee % boundary - face interior only 
    cc % corners and edges - ignored  

    %% faces 
    %    o---4---o
    %    |       |                       1,2 --> x=0,1
    %    1       2    face numbers ...   3,4 --> y=0,1
    %    |       |                       5,6 --> z=0,1
    %    o---3---o
    f1, f1s
    f2, f2s
    f3, f3s
    f4, f4s
    f5, f5s
    f6, f6s


    % Prolongation 
    p_h_1d
    r_h_1d
    r_hby2
    Ph      % interpolation from this element to its 4/8 children
    Rh
    Pp      % interpolation from this element to its 2p version		
    Pint
    Rint
    pcf_idx  % child face indices on interpolated parent volume
  end % properties 
    
  methods
    function elem = refel(d, order)
      % Setup a d-dimensional reference element 
      % order = polynomial order of the reference element
      % disp('refel');
      
      elem.dim  = d;
      elem.N    = order;
      elem.p    = order + 1;
      
      [elem.D, r] = mgps.refel.cheb(elem.N);
      elem.r = r(end:(-1):1);

      elem.D2 = elem.D^2; 
      I  = eye(elem.p);

      % interpolation operators
      r_hby2      = [0.5*(elem.r(2:order) - 1); 0.5*(elem.r(2:order) + 1)];
      r_hby2      = r_hby2(end:(-1):1);

      r_hby2_chl  = [0.5*(elem.r(1:order) - 1); 0.5*(elem.r + 1)];
      r_hby2_chl  = r_hby2_chl(end:(-1):1);
      
      [~, r_2p]   = mgps.refel.cheb (2*elem.N);
      
      p = elem.p;
      p2 = elem.p*elem.p;
      p3 = elem.p*elem.p*elem.p;

      Vr     = zeros (order-1, order-1);
      % elem.Vg     = zeros (order+1, order+1);
            
      Vph     = zeros (order-1, 2*(order-1));
			      
      for i=1:(elem.p-2)
        Vr(i,:)     = mgps.basis.polynomial (elem.r(2:order), 0, 0, i-1);
      %           elem.gradVr(i,:) = homg.basis.gradient (elem.r, 0, 0, i-1);
                
      %           elem.Vg(i,:)     = homg.basis.polynomial (elem.g, 0, 0, i-1);
      %           elem.gradVg(i,:) = homg.basis.gradient (elem.g, 0, 0, i-1);
                
        Vph(i,:)    = mgps.basis.polynomial (r_hby2, 0, 0, i-1);
%         Vpp(i,:)    = mgps.basis.polynomial (r_2p,   0, 0, i-1);
% 
%         Vph_par(i,:)    = mgps.basis.polynomial (elem.r, 0, 0, i-1);
%         Vph_chl(i,:)    = mgps.basis.polynomial (r_hby2_chl, 0, 0, i-1);
      end
        
      %       elem.Dr     = transpose(elem.Vr \ elem.gradVr);
            
      %       elem.Dg     = transpose(elem.Vr \ elem.gradVg);
            
      %       iVr         = elem.Vr \ eye(order+1);
			% 			iVg         = elem.Vg \ eye(order+1);
            
      %       elem.q1d         = transpose (elem.Vr \ elem.Vg);  
            
			p_h_1d      = transpose (Vr \ Vph); % flipud(p_h_1d);
      r_h_1d      = transpose (Vph \ Vr); % r_h_1d = flipud(r_h_1d);

%       p_h_1d_int      = transpose (Vph_par \ Vph_chl); flipud(p_h_1d_int);
%       r_h_1d_int      = transpose (Vph_chl \ Vph_par); flipud(r_h_1d_int);
%       p_p_1d      = transpose (Vr \ Vpp);  flipud(p_p_1d);
      
      elem.p_h_1d = p_h_1d;
      elem.r_hby2 = r_hby2;
      elem.r_h_1d = r_h_1d; 

      if (d == 2)
        elem.Dx = kron(I, elem.D);
        elem.Dy = kron(elem.D, I);mm

        elem.n = p2;

        elem.L = kron(I,elem.D2) + kron(elem.D2,I);
      else 
        % indices ...
        [i,j,k] = ndgrid(2:p-1, 2:p-1, 2:p-1);
        idx = sub2ind([p, p, p], i, j, k);
        elem.ii = sort( idx(:) );
        %%  faces 
        % f1 -> x=0
        [i,j,k] = ndgrid(1, 2:p-1, 2:p-1);
        idx = sub2ind([p, p, p], i, j, k);
        elem.f1 =  sort( idx(:) );
        % f2 -> x=1
        [i,j,k] = ndgrid(p, 2:p-1, 2:p-1);
        idx = sub2ind([p, p, p], i, j, k);
        elem.f2 =  sort( idx(:) );
        % f3 -> y=0
        [i,j,k] = ndgrid(2:p-1, 1, 2:p-1);
        idx = sub2ind([p, p, p], i, j, k);
        elem.f3 =  sort( idx(:) );        
        % f4 -> y=1
        [i,j,k] = ndgrid(2:p-1, p, 2:p-1);
        idx = sub2ind([p, p, p], i, j, k);
        elem.f4 =  sort( idx(:) );                
        % f5 --> z=0
        [i,j,k] = ndgrid(2:p-1, 2:p-1, 1);
        idx = sub2ind([p, p, p], i, j, k);
        elem.f5 =  sort( idx(:) );                
        % f6 --> z=1
        [i,j,k] = ndgrid(2:p-1, 2:p-1, p);
        idx = sub2ind([p, p, p], i, j, k);
        elem.f6 =  sort( idx(:) ); 

        elem.nnf = length(elem.f6);

        % face ordering - f1 f3 f2 f5 f4 f6
        elem.f1s = 1:elem.nnf;
        elem.f3s = (elem.nnf+1):(2*elem.nnf);
        elem.f2s = (2*elem.nnf+1):(3*elem.nnf);
        elem.f5s = (3*elem.nnf+1):(4*elem.nnf);
        elem.f4s = (4*elem.nnf+1):(5*elem.nnf);
        elem.f6s = (5*elem.nnf+1):(6*elem.nnf);

        % face ordering - f1 f3 f2 f5 f4 f6 
        elem.ee = [elem.f1; elem.f3; elem.f2; elem.f5; elem.f4; elem.f6];
        elem.cc = setdiff(1:p3, [elem.ii; elem.ee])';

        elem.n = p3;
        elem.ni = length(elem.ii);
        elem.ne = length(elem.ee);

        % elem.Ph = kron(p_h_1d, p_h_1d);

        elem.I   = kron(kron(I, I), I);
        elem.Dx  = kron(kron(I, I), elem.D);
        elem.Dy  = kron(kron(I, elem.D), I);
        elem.Dz  = kron(kron(elem.D, I), I);
        elem.Dxx = kron(kron(I, I), elem.D2); 
        elem.Dxy = elem.Dx * elem.Dy;
        elem.Dxz = elem.Dx * elem.Dz;
        elem.Dyy = kron(kron(I, elem.D2), I);
        elem.Dyz = elem.Dy * elem.Dz;
        elem.Dzz = kron(kron(elem.D2, I), I);
        elem.L = elem.Dxx + elem.Dyy + elem.Dzz;

        % face interpolations ...
        p1 = p_h_1d(1:end/2,:);
        p2 = p_h_1d(end/2+1:end,:);
        r1 = r_h_1d(:,1:end/2);
        r2 = r_h_1d(:,end/2+1:end);
        elem.Ph{1} = kron(p1, p1);
        elem.Ph{3} = kron(p2, p1);
        elem.Ph{2} = kron(p1, p2);
        elem.Ph{4} = kron(p2, p2);
        elem.Rh{1} = kron(r1, r1);
        elem.Rh{3} = kron(r2, r1);
        elem.Rh{2} = kron(r1, r2);
        elem.Rh{4} = kron(r2, r2);
%         elem.Pp = kron(p_p_1d, p_p_1d);
        %% for prolongation 
%         elem.Pint = kron(kron(p_h_1d_int,p_h_1d_int),p_h_1d_int);
        %------
        % dP = decomposition(elem.Pint);
        % elem.Rint = @(u) dP \ u;
%         elem.Rint = transpose(elem.Pint); %kron(kron(r_h_1d_int,r_h_1d_int),r_h_1d_int);%pinv( elem.Pint);
        %------
        elem.compute_face_nodes();       
      end 
    end % refel constructor

    function d = dxx(e, jac)
      d = e.Dxx/jac/jac;
    end 

    function d = dyy(e, jac)
      d = e.Dyy/jac/jac;
    end 

    function d = dzz(e, jac)
      d = e.Dzz/jac/jac;
    end
    
    function d = dxy(e, jac)
      d = e.Dxy/jac/jac;
    end 

    function d = dyz(e, jac)
      d = e.Dyz/jac/jac;
    end 

    function d = dxz(e, jac)
      d = e.Dxz/jac/jac;
    end

    function d = dx(e, jac)
      d = e.Dx/jac;
    end

    function d = dy(e, jac)
      d = e.Dy/jac;
    end

    function d = dz(e, jac)
      d = e.Dz/jac;
    end

    function d = b(e, jac)
      d = e.I;
    end

    function idx = compute_face_nodes(r)
      % gets child face indices for parent-level volume.
      pp = 2*r.p - 1;
      %=========== Child 1 =========== 
      [i,j,k] = ndgrid(1, 2:r.p-1, 2:r.p-1); % 1,1
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{1,1} = sort(idx(:));
      [i,j,k] = ndgrid(r.p, 2:r.p-1, 2:r.p-1); % 1,2
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{1,2} = sort(idx(:));
      [i,j,k] = ndgrid(2:r.p-1, 1, 2:r.p-1); % 1,3
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{1,3} = sort(idx(:));
      [i,j,k] = ndgrid(2:r.p-1, r.p, 2:r.p-1); % 1,4
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{1,4} = sort(idx(:));
      [i,j,k] = ndgrid(2:r.p-1, 2:r.p-1, 1); % 1,5
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{1,5} = sort(idx(:));
      [i,j,k] = ndgrid(2:r.p-1, 2:r.p-1, r.p); % 1,6
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{1,6} = sort(idx(:));
      %=========== Child 2 ===========
      [i,j,k] = ndgrid(r.p, 2:r.p-1, 2:r.p-1); % 2,1
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{2,1} = sort(idx(:));
      [i,j,k] = ndgrid(pp, 2:r.p-1, 2:r.p-1); % 2,2
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{2,2} = sort(idx(:));
      [i,j,k] = ndgrid((r.p+1):(pp-1), 1, 2:r.p-1); % 2,3
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{2,3} = sort(idx(:));
      [i,j,k] = ndgrid((r.p+1):(pp-1), r.p, 2:r.p-1); % 2,4
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{2,4} = sort(idx(:));
      [i,j,k] = ndgrid((r.p+1):(pp-1), 2:r.p-1, 1); % 2,5
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{2,5} = sort(idx(:));
      [i,j,k] = ndgrid((r.p+1):(pp-1), 2:r.p-1, r.p); % 2,6
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{2,6} = sort(idx(:));
      % Child 3  - (r.p+1):(pp-1)
      [i,j,k] = ndgrid(1, (r.p+1):(pp-1), 2:r.p-1); % 3,1
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{3,1} = sort(idx(:));
      [i,j,k] = ndgrid(r.p,(r.p+1):(pp-1), 2:r.p-1); % 3,2
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{3,2} = sort(idx(:));
      [i,j,k] = ndgrid(2:r.p-1, r.p, 2:r.p-1); % 3,3
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{3,3} = sort(idx(:));
      [i,j,k] = ndgrid(2:r.p-1, pp, 2:r.p-1); % 3,4
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{3,4} = sort(idx(:));
      [i,j,k] = ndgrid(2:r.p-1, (r.p+1):(pp-1), 1); % 3,5
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{3,5} = sort(idx(:));
      [i,j,k] = ndgrid(2:r.p-1, (r.p+1):(pp-1), r.p); % 3,6
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{3,6} = sort(idx(:));
      % Child 4
      [i,j,k] = ndgrid(r.p, (r.p+1):(pp-1), 2:r.p-1); % 4,1
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{4,1} = sort(idx(:));
      [i,j,k] = ndgrid(pp, (r.p+1):(pp-1), 2:r.p-1); % 4,2
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{4,2} = sort(idx(:));
      [i,j,k] = ndgrid((r.p+1):(pp-1), r.p, 2:r.p-1); % 4,3
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{4,3} = sort(idx(:));
      [i,j,k] = ndgrid((r.p+1):(pp-1), pp, 2:r.p-1); % 4,4
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{4,4} = sort(idx(:));
      [i,j,k] = ndgrid((r.p+1):(pp-1), (r.p+1):(pp-1), 1); % 4,5
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{4,5} = sort(idx(:));
      [i,j,k] = ndgrid((r.p+1):(pp-1), (r.p+1):(pp-1), r.p); % 4,6
      idx = sub2ind([pp, pp, pp], i, j, k);
      r.pcf_idx{4,6} = sort(idx(:));
          %=========== Child 5 =========== 
          [i,j,k] = ndgrid(1, 2:r.p-1, (r.p+1):(pp-1)); % 1,1
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{5,1} = sort(idx(:));
          [i,j,k] = ndgrid(r.p, 2:r.p-1, (r.p+1):(pp-1)); % 1,2
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{5,2} = sort(idx(:));
          [i,j,k] = ndgrid(2:r.p-1, 1, (r.p+1):(pp-1)); % 1,3
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{5,3} = sort(idx(:));
          [i,j,k] = ndgrid(2:r.p-1, r.p, (r.p+1):(pp-1)); % 1,4
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{5,4} = sort(idx(:));
          [i,j,k] = ndgrid(2:r.p-1, 2:r.p-1, r.p); % 1,5
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{5,5} = sort(idx(:));
          [i,j,k] = ndgrid(2:r.p-1, 2:r.p-1, pp); % 1,6
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{5,6} = sort(idx(:));
          %=========== Child 6 ===========
          [i,j,k] = ndgrid(r.p, 2:r.p-1, (r.p+1):(pp-1)); % 2,1
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{6,1} = sort(idx(:));
          [i,j,k] = ndgrid(pp, 2:r.p-1, (r.p+1):(pp-1)); % 2,2
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{6,2} = sort(idx(:));
          [i,j,k] = ndgrid((r.p+1):(pp-1), 1, (r.p+1):(pp-1)); % 2,3
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{6,3} = sort(idx(:));
          [i,j,k] = ndgrid((r.p+1):(pp-1), r.p, (r.p+1):(pp-1)); % 2,4
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{6,4} = sort(idx(:));
          [i,j,k] = ndgrid((r.p+1):(pp-1), 2:r.p-1, r.p); % 2,5
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{6,5} = sort(idx(:));
          [i,j,k] = ndgrid((r.p+1):(pp-1), 2:r.p-1, pp); % 2,6
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{6,6} = sort(idx(:));
          %========== Child 7 ================ 
          [i,j,k] = ndgrid(1, (r.p+1):(pp-1), (r.p+1):(pp-1)); % 3,1
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{7,1} = sort(idx(:));
          [i,j,k] = ndgrid(r.p,(r.p+1):(pp-1), (r.p+1):(pp-1)); % 3,2
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{7,2} = sort(idx(:));
          [i,j,k] = ndgrid(2:r.p-1, r.p, (r.p+1):(pp-1)); % 3,3
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{7,3} = sort(idx(:));
          [i,j,k] = ndgrid(2:r.p-1, pp, (r.p+1):(pp-1)); % 3,4
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{7,4} = sort(idx(:));
          [i,j,k] = ndgrid(2:r.p-1, (r.p+1):(pp-1), r.p); % 3,5
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{7,5} = sort(idx(:));
          [i,j,k] = ndgrid(2:r.p-1, (r.p+1):(pp-1), pp); % 3,6
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{7,6} = sort(idx(:));
          % Child 8
          [i,j,k] = ndgrid(r.p, (r.p+1):(pp-1), (r.p+1):(pp-1)); % 4,1
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{8,1} = sort(idx(:));
          [i,j,k] = ndgrid(pp, (r.p+1):(pp-1), (r.p+1):(pp-1)); % 4,2
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{8,2} = sort(idx(:));
          [i,j,k] = ndgrid((r.p+1):(pp-1), r.p, (r.p+1):(pp-1)); % 4,3
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{8,3} = sort(idx(:));
          [i,j,k] = ndgrid((r.p+1):(pp-1), pp, (r.p+1):(pp-1)); % 4,4
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{8,4} = sort(idx(:));
          [i,j,k] = ndgrid((r.p+1):(pp-1), (r.p+1):(pp-1), r.p); % 4,5
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{8,5} = sort(idx(:));
          [i,j,k] = ndgrid((r.p+1):(pp-1), (r.p+1):(pp-1), pp); % 4,6
          idx = sub2ind([pp, pp, pp], i, j, k);
          r.pcf_idx{8,6} = sort(idx(:));
      end

  end % methods

  methods (Static)
    function [D, x] = cheb(N)
      if (N==0) 
        D=0; 
        x=1; 
        return 
      end

      x = cos(pi*(0:N)/N)'; 
      c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
      X = repmat(x,1,N+1);
      dX = X-X';                  
      D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
      D  = D - diag(sum(D'));                 % diagonal entries
    end
  end % private methods
    
end

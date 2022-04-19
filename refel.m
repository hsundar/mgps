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

    % interpolation 
    q1d
    
    G2C

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
    Ph      % interpolation from this element to its 4/8 children
    Pp      % interpolation from this element to its 2p version		
  end % properties 
    
  methods
    function elem = refel(d, order)
      % Setup a d-dimensional reference element 
      % order = polynomial order of the reference element
      
      elem.dim  = d;
      elem.N    = order;
      elem.p    = order + 1;
      
      [elem.D, r] = refel.cheb(elem.N);
      elem.r = r(end:(-1):1);

      elem.D2 = elem.D^2; 
      I  = eye(elem.p);

      % interpolation operators
      % r_hby2      = [0.5*(elem.r - 1); 0.5*(elem.r(2:end) + 1)];
      % [_, r_2p]   = refel.cheb (2*elem.N);
      Nrp = N+1;

      Vr     = zeros (order+1, order+1);
      Vg     = zeros (order+1, order+1);
      
      for i=1:Nrp
        Vr(i,:)           = mgps.basis.polynomial (elem.r, 0, 0, i-1);
        Vg(i,:)           = mgps.basis.polynomial (elem.g, 0, 0, i-1);
        
        %Vph(i,:)          = homg.basis.polynomial (r_hby2, 0, 0, i-1);
        %Vpp(i,:)          = homg.basis.polynomial (r_2p,   0, 0, i-1);
      end
    
      elem.q1d         = transpose (Vr \ Vg);  

      
      p = elem.p;
      p2 = elem.p*elem.p;
      p3 = elem.p*elem.p*elem.p;
      
      if (d == 2)
        elem.Dx = kron(I, elem.D);
        elem.Dy = kron(elem.D, I);

        elem.n = p2;

        elem.L = kron(I,elem.D2) + kron(elem.D2,I);

        elem.Q  = kron(elem.q1d, elem.q1d) ;
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
        % g2c 
        Qf  = kron(elem.q1d, elem.q1d);
        elem.G2C = sparse(elem.ne, 6*elem.nnf);

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


classdef grid < handle
    %GRID A single grid in a multigrid heirarchy
    
    properties
        level
        is_finest
        
        order
        Mesh
        Coarse  % handle to coarse grid
        debug
        dbg_spaces
        
        Boundary
        BoundaryValues
        BoundaryPoints
    end % properties
    
    methods
        function grid = grid(mesh, order, coarse)
            if ((nargin < 3) || isempty(coarse))
                grid.level = 0;
                grid.Coarse = [];
            else
                grid.level = coarse.level + 1;
                grid.Coarse = coarse;
            end
            
            grid.dbg_spaces = '      ';
            grid.dbg_spaces = grid.dbg_spaces(1:end-4*grid.level);
            grid.debug = 0;
            
            grid.Mesh = mesh;
            grid.order = order;
            
            % mesh.set_order(order);
%             if (~ isempty(grid.Coarse) )
%                 grid.P = grid.Coarse.Mesh.assemble_interpolation(order);
%                 grid.R = grid.P';
%             end
            bdy_idx = [];
            num_elems  = prod(grid.Mesh.nelems);
            r = mgps.refel(3,order);
            nnf = r.nnf;
            bdy_pts = [];
            for e=1:num_elems
                pts = grid.Mesh.element_nodes(e, r);
                fid = mesh.get_global_faces(e);
                for f=1:length(fid)
                    bdy_idx = [bdy_idx, ((fid(f)-1)*nnf+1):(fid(f)*nnf)];
                end
                bdy_pts = [bdy_pts; pts(r.f1,:); pts(r.f3,:); pts(r.f2,:); pts(r.f5,:); pts(r.f4,:); pts(r.f6,:)];
            end
            % for Dirichlet boundary conditions
            grid.Boundary = bdy_idx;  % mesh.get_boundary_node_indices(order);
            grid.BoundaryPoints = bdy_pts;
            % grid.BoundaryValues = feval(gx, bdy_pts(:,1), pts(:,2), pts(:,3));

            grid.is_finest       = false;            
        end
        
        function assemble_operators(grid, op, mu, rhs, bdy_fx)
            % assemble for this level ...
            grid.Mesh.initialize_leaves(grid.order, op, rhs, mu);
            
            grid.BoundaryValues = feval(bdy_fx, grid.BoundaryPoints(:,1), grid.BoundaryPoints(:,2), grid.BoundaryPoints(:,3));
            
            % propagate to lower grids
            if (~ isempty(grid.Coarse) )
                if isnumeric(mu)
                    harmonic = 0;
                    if (grid.Mesh.dim == 2)
                        mu2 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2));
                        if (harmonic)
                            mu_coarse = 4 ./ ( 1./mu2(1:2:end, 1:2:end) + 1./mu2(2:2:end, 1:2:end) + 1./mu2(1:2:end, 2:2:end) + 1./mu2(2:2:end, 2:2:end) );
                        else
                            mu_coarse = 0.25*(mu2(1:2:end, 1:2:end) + mu2(2:2:end, 1:2:end) + mu2(1:2:end, 2:2:end) + mu2(2:2:end, 2:2:end));
                        end
                    else
                        mu3 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2), grid.Mesh.nelems(3));
                        if (harmonic)
                            mu_coarse = 8 ./ ( 1./mu3(1:2:end, 1:2:end, 1:2:end) + 1./mu3(2:2:end, 1:2:end, 1:2:end) ...
                                + 1./mu3(1:2:end, 2:2:end, 1:2:end) + 1./mu3(2:2:end, 2:2:end, 1:2:end) ...
                                + 1./mu3(1:2:end, 1:2:end, 2:2:end) + 1./mu3(2:2:end, 1:2:end, 2:2:end) ...
                                + 1./mu3(1:2:end, 2:2:end, 2:2:end) + 1./mu3(2:2:end, 2:2:end, 2:2:end) );
                        else
                            mu_coarse = 0.125*(mu3(1:2:end, 1:2:end, 1:2:end) + mu3(2:2:end, 1:2:end, 1:2:end) + mu3(1:2:end, 2:2:end, 1:2:end) + mu3(2:2:end, 2:2:end, 1:2:end) + ...
                                mu3(1:2:end, 1:2:end, 2:2:end) + mu3(2:2:end, 1:2:end, 2:2:end) + mu3(1:2:end, 2:2:end, 2:2:end) + mu3(2:2:end, 2:2:end, 2:2:end) );
                        end
                    end
                    grid.Coarse.assemble_operators(op, mu_coarse, rhs, bdy_fx);
                else
                    grid.Coarse.assemble_operators(op, mu, rhs, bdy_fx);
                end
            end
        end
                
        % compute the residual
        function r = residual(grid, rhs, u)
            % function r = residual(grid, u, rhs)
            if ( nargin < 2 )
                rhs = grid.L;
            end
            if ( nargin < 3 )
                u = zeros(size(rhs));
            end
            
            r = rhs - grid.K*u;
        end
                
        function [u, rr, iter] = solve_pcg(grid, num_vcyc, smoother, v1, v2, rhs, u)
            % disp('setting smoother');
            grid.set_smoother(smoother);
            
            % disp('computing initial residual');
            r = grid.residual(rhs, u);
            rho = zeros(size(u));
            % disp('outer v-cycle');
            rho = grid.vcycle(v1, v2, r, rho);
            p = rho;
            disp(['Initial residual is ' num2str(norm(r))]);
            disp('------------------------------------------');
            r0 = norm(r);
            for i=1:num_vcyc
                % disp(['inner v-cycle: ' num2str(i)]);
                h = grid.K * p;
                rho_res = dot (rho, r);
                alpha = rho_res / dot ( p, h );
                u = u + alpha*p;
                r = r - alpha*h;
                
                % rho_res_prev = rho_res;
                
                disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e')]);
                if (norm(r)/r0 < 1e-8)
                    iter = i;
                    rr = norm(r)/r0;
                    return;
                end
                
                % precondition ..
                rho = zeros(size(u)); % needed ?
                rho = grid.vcycle(v1, v2, r, rho);
                
                beta = dot(rho, r) / rho_res ;
                p = rho + beta*p;
            end
            disp('------------------------------------------');
            iter = num_vcyc;
            rr = norm(r)/r0;
        end
          
        function [u, rr, iter] = solve(grid, num_vcyc, smoother, v1, v2, rhs, u)
            grid.set_smoother(smoother);
            
            r = grid.residual(rhs, u);
            
            disp(['Initial residual is ' num2str(norm(r))]);
            disp('------------------------------------------');
            r0 = norm(r);
            
            for i=1:num_vcyc
                u = grid.vcycle(v1, v2, rhs, u);
                r = grid.residual(rhs, u);
                disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e')]);
                if (norm(r)/r0 < 1e-8)
                    iter = i;
                    rr = norm(r)/r0;
                    return;
                end
            end
            disp('------------------------------------------');
            iter = num_vcyc;
            rr = norm(r)/r0;
        end
        
        function u = vcycle(grid, v1, v2, rhs, u)
            % function u = vcycle(grid, v1, v2, rhs, u)
            % solve system using initial guess u, given rhs
            % with v1 pre and v2 post-smoothing steps
            % disp(['CG vcycle: order ' num2str(grid.Mesh.order) ', nelems: ' num2str(grid.Mesh.nelems(1)) 'X' num2str(grid.Mesh.nelems(2))]);
            
            if ( isempty( grid.Coarse ) )
                if (grid.linear_smoother)
                    u = grid.K_lin \ rhs;
                else
                    u = grid.K \ rhs;
%                    u = grid.K \ (grid.M * rhs);
                end
                
                return;
            end
            
            % 1. pre-smooth
            u = grid.smooth ( v1, rhs, u );
                
            % 2. compute residual
            if (grid.linear_smoother && ~grid.is_finest)
                % disp('linear residual');
                res = grid.residual_lin(rhs, u);
            else
                % disp('high-order residual');
                res = grid.residual(rhs, u);
            end
                
           % 3. restrict
           res_coarse = grid.R * res;
           res_coarse(grid.Coarse.Boundary) = 0;
            
           % 4. ---------- recurse -----------
           u_corr_coarse = grid.Coarse.vcycle(v1, v2, res_coarse, zeros(size(res_coarse)));
            
           % 5. prolong and correct
           u = u + grid.P * u_corr_coarse;
           
           % 6. post-smooth
           u = grid.smooth ( v2, rhs, u );
           % grid.plot_spectrum(u, 'g', rhs);
            
       end % v-cycle
       
        
        function set_coeff(grid, mu)
            grid.Mesh.set_coeff (mu) ;
            if (~ isempty(grid.Coarse) )
                grid.Coarse.Mesh.set_coeff (mu) ;
            end
        end
        
                
        function u0 = get_u0(grid)
            if (grid.debug)
                if ( isempty( grid.k_evec ) )
                    [grid.k_evec, ~] = eigs(grid.K, grid.M, 80, 'BE');
                end
                n = size(grid.k_evec, 2);
                lam = ones(n,1);
                % lam(1:n/4) = 1;
                u0 = grid.k_evec*lam;
            else
                u0 = rand(size(grid.L()));
                u0(grid.Boundary) = 0;
            end
        end
                
    end %methods
    
end %classdef


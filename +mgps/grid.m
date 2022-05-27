classdef grid < handle
    %GRID A single grid in a multigrid heirarchy
    
    properties
        level
        is_finest
        
        order
        Mesh
        Coarse  % handle to coarse grid
        Fine    % handle to fine Grid
        debug
        dbg_spaces
        
        smoother 
        jacobi_omega
        jacobi_invdiag

        Boundary
        BoundaryValues
        BoundaryPoints

        mergemap

        pfac
    end % properties
    
    methods
        function grid = grid(mesh, order, coarse)
            if ((nargin < 3) || isempty(coarse))
                grid.level = 0;
                grid.Coarse = [];
            else
                grid.level = coarse.level + 1;
                grid.Coarse = coarse;
                coarse.Fine = grid;
            end
            
            grid.dbg_spaces = '      ';
            grid.dbg_spaces = grid.dbg_spaces(1:end-4*grid.level);
            grid.debug = 0;
            grid.pfac = 1;
            
            grid.Mesh = mesh;
            grid.order = order;
            grid.smoother = 'jacobi';
            grid.jacobi_omega = 6/7;
            
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

            grid.mergemap        = mgps.mergemaps;
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
        function r = residual(grid, u)
            % function r = residual(grid, u, rhs)
%             if ( nargin < 2 )
%                 u = grid.get_u0();
%             end
%             
            % r = rhs - grid.K*u;
            r = grid.Mesh.trace_residual(u);

            % grid.plot_skel(r); keyboard
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
          
        function [u, rr, iter] = solve(grid, num_vcyc, smoother, v1, v2, u)
            grid.set_smoother(smoother);
            
            r = grid.residual(u);
            
            disp(['Initial residual is ' num2str(norm(r))]);
            disp('------------------------------------------');
            r0 = norm(r);
            
            for i=1:num_vcyc
%                 plot(u);
%                 pause
                u = grid.vcycle(v1, v2, u);
                r = grid.residual(u);
%                 grid.plot_skel(r);
%                 pause
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
        
        function u = vcycle(grid, v1, v2, u)
            % function u = vcycle(grid, v1, v2, rhs, u)
            % solve system using initial guess u, given rhs
            % with v1 pre and v2 post-smoothing steps
            % disp(['CG vcycle: order ' num2str(grid.Mesh.order) ', nelems: ' num2str(grid.Mesh.nelems(1)) 'X' num2str(grid.Mesh.nelems(2))]);
            
            if ( isempty( grid.Coarse ) )
                % u = grid.K \ rhs;
                %tstart = tic;
                u = grid.smooth(20, u);
                %tend = toc(tstart);
                %fprintf('Coarse(st) time: %f \n', tend);
                return;
            end
            
            % 1. pre-smooth
%             if grid.is_finest
%             tstart = tic;
%             end
            u = grid.smooth ( v1, u );
%             if grid.is_finest
%             tend = toc(tstart);
%             fprintf('Smooth time: %f \n', 2*tend);
%             
%             tstart = tic;
%             end
            % 2. compute residual
            res = grid.residual(u);
                
            % 3. restrict
            %res_coarse = 
            grid.restrict(res);

%             if grid.is_finest
%             tend = toc(tstart);
%             fprintf('Restrict time: %f \n', tend);
%             end
            % res_coarse = grid.R * res;
            % res_coarse(grid.Coarse.Boundary) = 0;
            
%             res_coarse = grid.Coarse.residual();
%             grid.Coarse.plot_skel(res_coarse);
%             pause;
            
            % 4. ---------- recurse -----------
%             if grid.is_finest
%             tstart = tic;
%             end
            u_corr_coarse = grid.Coarse.vcycle(v1, v2, grid.Coarse.get_u0());
%             if grid.is_finest
%             tend = toc(tstart);
%             fprintf('Coarse time: %f \n', tend);
%             end

            % 5. prolong and correct
            
%             if grid.is_finest
%             tstart = tic;
%             end
            uc = grid.prolong( u_corr_coarse );
            u = u - grid.pfac*uc; 
%             if grid.is_finest
%             tend = toc(tstart);
%             fprintf('Prolong time: %f \n', tend);
%             end

            % 6. post-smooth
            u = grid.smooth ( v2, u );
        
%             grid.plot_skel(u);
%             pause;
        end % v-cycle
       
       % smoothers
       function u = smooth (grid, v, u)
        switch(grid.smoother)
            case 'jacobi',
                u = grid.smoother_jacobi(v, u);
                return;
            case 'l1_jac',
                u = grid.smoother_l1_jacobi(v, u);
                return;
            case 'blk_jac',
                u = grid.smoother_block_jacobi(v, u);
                return;
            case 'gs',
                grid.sor_omega = 1.0;
                u = grid.smoother_gauss_seidel(v, u);
                return;
            case 'chebssor'
                u = grid.smoother_chebyshev_ssor (v, u);
                return;
            case 'chebyshev2'
                u = grid.smoother_chebyshev2 (v, u);
                return;
            case 'chebyshev',
                u = grid.smoother_chebyshev(v, u);
                return;
            case 'sor',
                u = grid.smoother_sor(v, u);
                return;
            case 'ssor',
                u = grid.smoother_sym_sor(v, u);
                return;
            case '2sr',
                u = grid.smoother_2sr(v, u);
                return;
            case 'hybrid',
                u = grid.smoother_hybrid(v, u);
                return;
            otherwise
                disp('ERROR: Unrecognized smoother type');
                return;
            end
        end
    
        function set_smoother(grid, sm)
            grid.smoother = sm;
            if (~ isempty(grid.Coarse) )
                grid.Coarse.set_smoother(sm);
            end
        end

        function u = smoother_jacobi (grid, v,  u)
            % standard jacobi smoother
            if ( isempty(grid.jacobi_invdiag) )
                D = grid.Mesh.trace_diagonal();
                grid.jacobi_invdiag = 1./D;
            end
            for i=1:v
              res  = grid.residual(u);  
              % grid.plot_skel(res);
            %   pause
              %~~~~ DEBUG ~~~~
              % plot(res);
%               nnf = grid.Mesh.refel.nnf;
%               fr = sqrt(nnf);
%               I = zeros(32,32);
%               for k=1:8
%                 for j=1:8
%                   I( ((k-1)*fr+1):(k*fr), ((j-1)*fr+1):(j*fr)) = reshape(res( ((8*(k-1) + j-1)*nnf+1):(8*(k-1) + j)*nnf),fr,fr) ; 
%                 end
%               end
%               imagesc(I); colorbar; axis equal;
%               pause
              %%~~~~~~~~~~~~~~  
              u = u - grid.jacobi_omega.*grid.jacobi_invdiag .*res;
%               r = norm(res);
%               disp([grid.dbg_spaces num2str(r)]);
              %norm(r)
            end
        end % jacobi

        %% add other smoothers 


        function set_coeff(grid, mu)
            grid.Mesh.set_coeff (mu) ;
            if (~ isempty(grid.Coarse) )
                grid.Coarse.Mesh.set_coeff (mu) ;
            end
        end
        
                
        function u0 = get_u0(grid)
            nnf = grid.Mesh.refel.nnf;
            num_elems  = prod(grid.Mesh.nelems);
            nelems = grid.Mesh.nelems;
            num_bdy = nnf * ( 3*num_elems + nelems(1)*nelems(2) + nelems(2)*nelems(3) + nelems(1)*nelems(3));
            u0 = zeros(num_bdy,1);
        end
                
        function r = prolong(grid, rc)
            num_elem_c = grid.Coarse.Mesh.nelems;
            num_elem_f = grid.Mesh.nelems;
            nnf = grid.Mesh.refel.nnf;
            forder = [1 3 2 5 4 6];
            % for all elems - coarse
            r = grid.get_u0();
            for k=1:num_elem_c(3)
                for j=1:num_elem_c(2)
                    for i=1:num_elem_c(1)
                        ep = sub2ind (num_elem_c, i, j, k);
                        fid_p = grid.Coarse.Mesh.get_global_faces(ep);
                        % use solution operator to get interior faces
                        bdy_idx = [];
                        % copy trace to element boundaries 
                        for f=1:length(fid_p)
                           bdy_idx = [bdy_idx, ((fid_p(f)-1)*nnf+1):(fid_p(f)*nnf)];
                        end
                        up = grid.Coarse.Mesh.S{ep}*[rc(bdy_idx); 1];
                        upf = grid.Mesh.refel.Pint * up;
                        %------------------------------ 
                        ce(1) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k-1);
                        ce(2) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k-1);
                        ce(3) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k-1);
                        ce(4) = sub2ind (num_elem_f, 2*i, 2*j, 2*k-1);
                        ce(5) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k);
                        ce(6) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k);
                        ce(7) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k);
                        ce(8) = sub2ind (num_elem_f, 2*i, 2*j, 2*k);
                        for c=1:8
                            fid_c = grid.Mesh.get_global_faces(ce(c));
                            for f=1:3 
                                % 3 shared faces with parent
                                ff = forder(grid.mergemap.ploc(c,f));
                                idx = grid.Mesh.refel.pcf_idx{c, grid.mergemap.ploc(c,f)};
                                r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf)) = r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf)) + upf(idx);
                                % 3 shared faces with siblings
                                ff = forder(grid.mergemap.sloc(c,f));
                                idx = grid.Mesh.refel.pcf_idx{c,grid.mergemap.sloc(c,f)};
                                r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf)) = r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf)) + upf(idx); 
                            end
                        end
                        
                    end % i
                end % j
            end % k
            % boundaries ?
        end % function prolong 

        function restrict(grid, r)
            num_elem_c = grid.Coarse.Mesh.nelems;
            num_elem_f = grid.Mesh.nelems;
            nnf = grid.Mesh.refel.nnf;
            nr = grid.Mesh.refel.p;
            forder = [1 3 2 5 4 6];

            % for all elems - coarse
            for k=1:num_elem_c(3)
                for j=1:num_elem_c(2)
                    for i=1:num_elem_c(1)
                        % rc = zeros(nnf*6,1);
                        %rcpp = zeros((2*grid.Mesh.refel.p - 1)^3,1);
                        ep = sub2ind (num_elem_c, i, j, k);

                        ce(1) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k-1);
                        ce(2) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k-1);
                        ce(3) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k-1);
                        ce(4) = sub2ind (num_elem_f, 2*i, 2*j, 2*k-1);
                        ce(5) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k);
                        ce(6) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k);
                        ce(7) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k);
                        ce(8) = sub2ind (num_elem_f, 2*i, 2*j, 2*k);
                        for c=1:8
                            fid_c = grid.Mesh.get_global_faces(ce(c));
                            bdy_idx = [];
                            for f=1:length(fid_c)
                                bdy_idx = [bdy_idx, ((fid_c(f)-1)*nnf+1):(fid_c(f)*nnf)];
                            end
                            rc = grid.Coarse.Mesh.S{ep}*[r(bdy_idx); 1];
                       
                            [ci, cj, ck] = ind2sub([2,2,2], c);
                            % x(((ci-1)*(nr-1)+1):ci*(nr-1)+1)
                            rcpp(((ci-1)*(nr-1)+1):ci*(nr-1)+1, ((cj-1)*(nr-1)+1):cj*(nr-1)+1, ((ck-1)*(nr-1)+1):ck*(nr-1)+1) = reshape(rc, [nr, nr, nr]);
                        end
                        rhs = grid.Mesh.refel.Rint * rcpp(:);
                        % rhs = transpose(grid.Mesh.refel.Pint) * rcpp(:);
                        grid.Coarse.Mesh.update_rhs(ep, rhs);
                    end % i
                end % j
            end % k
            
        end % function restrict 

        function restrict_fast(grid, r)
            num_elem_c = grid.Coarse.Mesh.nelems;
            num_elem_f = grid.Mesh.nelems;
            nnf = grid.Mesh.refel.nnf;
            forder = [1 3 2 5 4 6];
            % for all elems - coarse
           %  rc = grid.Coarse.get_u0();
            for k=1:num_elem_c(3)
                for j=1:num_elem_c(2)
                    for i=1:num_elem_c(1)
                        rc = zeros(nnf*6,1);
                        ep = sub2ind (num_elem_c, i, j, k);
                        % fid_p = grid.Coarse.Mesh.get_global_faces(ep);
                        % for f=1:length(fid)
                        %     bdy_idx_p = [bdy_idx_p, ((fid(f)-1)*nnf+1):(fid(f)*nnf)];
                        % end
                        ce(1) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k-1);
                        ce(2) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k-1);
                        ce(3) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k-1);
                        ce(4) = sub2ind (num_elem_f, 2*i, 2*j, 2*k-1);
                        ce(5) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k);
                        ce(6) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k);
                        ce(7) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k);
                        ce(8) = sub2ind (num_elem_f, 2*i, 2*j, 2*k);
                        for c=1:8
                            fid_c = grid.Mesh.get_global_faces(ce(c));
                            for f=1:3 % 3 shared faces with parent
                                % restrict/prolong
                                ff = forder(grid.mergemap.ploc(c,f));
                                rc(((ff-1)*nnf+1):(ff*nnf)) = rc(((ff-1)*nnf+1):(ff*nnf)) + transpose(grid.Mesh.refel.Ph{grid.mergemap.igrid(c,f)}) * r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf));
                            end
                        end
                        % grid.Mesh.refel.n
                        rhs = grid.Coarse.Mesh.S{ep} * [rc; 1];
                        grid.Coarse.Mesh.update_rhs(ep, rhs);
                    end % i
                end % j
            end % k
            
        end % function restrict_fast 

        function rc = restrict_old(grid, r)
            num_elem_c = grid.Coarse.Mesh.nelems;
            num_elem_f = grid.Mesh.nelems;
            nnf = grid.Mesh.refel.nnf;
            forder = [1 3 2 5 4 6];
            % for all elems - coarse
            rc = grid.Coarse.get_u0();
            for k=1:num_elem_c(3)
                for j=1:num_elem_c(2)
                    for i=1:num_elem_c(1)
                        ep = sub2ind (num_elem_c, i, j, k);
                        fid_p = grid.Coarse.Mesh.get_global_faces(ep);
                        % for f=1:length(fid)
                        %     bdy_idx_p = [bdy_idx_p, ((fid(f)-1)*nnf+1):(fid(f)*nnf)];
                        % end
                        ce(1) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k-1);
                        ce(2) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k-1);
                        ce(3) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k-1);
                        ce(4) = sub2ind (num_elem_f, 2*i, 2*j, 2*k-1);
                        ce(5) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k);
                        ce(6) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k);
                        ce(7) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k);
                        ce(8) = sub2ind (num_elem_f, 2*i, 2*j, 2*k);
                        for c=1:8
                            fid_c = grid.Mesh.get_global_faces(ce(c));
                            for f=1:3 % 3 shared faces with parent
                                % restrict/prolong
                                ff = forder(grid.mergemap.ploc(c,f));
                                rc(((fid_p(ff)-1)*nnf+1):(fid_p(ff)*nnf)) = rc(((fid_p(ff)-1)*nnf+1):(fid_p(ff)*nnf)) + transpose(grid.Mesh.refel.Ph{grid.mergemap.igrid(c,f)}) * r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf));
                            end
                        end
                    end % i
                end % j
            end % k
            %% zero out bdy ?
        end % function restrict_old 

        function r = prolong_old(grid, rc)
            num_elem_c = grid.Coarse.Mesh.nelems;
            num_elem_f = grid.Mesh.nelems;
            nnf = grid.Mesh.refel.nnf;
            forder = [1 3 2 5 4 6];
            % for all elems - coarse
            r = grid.get_u0();
            for k=1:num_elem_c(3)
                for j=1:num_elem_c(2)
                    for i=1:num_elem_c(1)
                        ep = sub2ind (num_elem_c, i, j, k);
                        fid_p = grid.Coarse.Mesh.get_global_faces(ep);
                        % use solution operator to get interior faces
                        % bdy_idx = [];
                        % % copy trace to element boundaries 
                        % for f=1:length(fid_p)
                        %   bdy_idx = [bdy_idx, ((fid_p(f)-1)*nnf+1):(fid_p(f)*nnf)];
                        % end
                        % up = grid.Coarse.Mesh.S{ep}*[rc(bdy_idx); 1];
                        % upf = grid.Mesh.refel.Pint * up;
                        %------------------------------ 
                        ce(1) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k-1);
                        ce(2) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k-1);
                        ce(3) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k-1);
                        ce(4) = sub2ind (num_elem_f, 2*i, 2*j, 2*k-1);
                        ce(5) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k);
                        ce(6) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k);
                        ce(7) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k);
                        ce(8) = sub2ind (num_elem_f, 2*i, 2*j, 2*k);
                        for c=1:8
                            fid_c = grid.Mesh.get_global_faces(ce(c));
                            for f=1:3 
                                % restrict/prolong
                                ff = forder(grid.mergemap.ploc(c,f));
                                r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf)) = grid.Mesh.refel.Ph{grid.mergemap.igrid(c,f)} * rc(((fid_p(ff)-1)*nnf+1):(fid_p(ff)*nnf));
                                % 3 shared faces with parent

                                % 3 shared faces with siblings
                            end
                        end
                        
                    end % i
                end % j
            end % k
        end % function prolong_old

        function plot_skel(grid, u)
          nelems = grid.Mesh.nelems;
          nnf = grid.Mesh.refel.nnf;
          fr = sqrt(nnf);

          % grid
          s = [nelems(1) nelems(2)]; % [y x]
          xrange = [-1 1]; % imagesc only needs the endpoints
          yrange = [-1 1];
          dx = diff(xrange)/(s(2)*fr -1);
          dy = diff(yrange)/(s(1)*fr -1);
          xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
          yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);

          %% x=const faces
          % x = 0
          subplot(3,3,1);
          I = zeros(nelems(2)*fr, nelems(3)*fr);
          for k=1:nelems(3)
            for j=1:nelems(2)
                idx = (k-1)*nelems(2)*(nelems(1)+1) + (j-1)*(nelems(1)+1) + 1;
                I( ((j-1)*fr+1):(j*fr),((k-1)*fr+1):(k*fr) ) = reshape(u((idx*nnf+1):(idx+1)*nnf),fr,fr) ; 
            end
          end
          imagesc(xrange,yrange,I); hold on; colorbar; axis equal; 
          h1 = mesh(xg,yg,zeros(s+1));
          h1.FaceColor = 'none';
          h1.EdgeColor = 'k';  

          title('x=1'); hold off
        %   offset = nelems(3)*nelems(2)*(nelems(1)/2) ;
          subplot(3,3,2);
          for k=1:nelems(3)
            for j=1:nelems(2)
                idx = (k-1)*nelems(2)*(nelems(1)+1) + (j-1)*(nelems(1)+1) + nelems(1)/2;
                I( ((j-1)*fr+1):(j*fr),((k-1)*fr+1):(k*fr) ) = reshape(u((idx*nnf+1):(idx+1)*nnf),fr,fr) ;  
            end
          end
          imagesc(xrange,yrange,I); colorbar; axis equal; hold on;
          h2 = mesh(xg,yg,zeros(s+1));
          h2.FaceColor = 'none';
          h2.EdgeColor = 'k';  
          title('x=mid'); hold off
        %   offset = nelems(3)*nelems(2)*(nelems(1)) ;
          subplot(3,3,3);
          for k=1:nelems(3)
            for j=1:nelems(2)
                idx = (k-1)*nelems(2)*(nelems(1)+1) + (j-1)*(nelems(1)+1) + nelems(1) - 1;
                I( ((j-1)*fr+1):(j*fr),((k-1)*fr+1):(k*fr) ) = reshape(u((idx*nnf+1):(idx+1)*nnf),fr,fr) ; 
            end
          end
          imagesc(xrange,yrange,I); colorbar; axis equal; hold on;
          h3 = mesh(xg,yg,zeros(s+1));
          h3.FaceColor = 'none';
          h3.EdgeColor = 'k';  
          title('x=end'); hold off
          offset = nelems(3)*nelems(2)*(nelems(1)+1) ;
          %% y=const faces
          I = zeros(nelems(1)*fr, nelems(3)*fr); 
          % y = 0
          subplot(3,3,4);
          for k=1:nelems(3)
            for j=1:nelems(2)
                idx = offset + (k-1)*(nelems(2)+1)*nelems(1)  + j-1;
                I( ((j-1)*fr+1):(j*fr),((k-1)*fr+1):(k*fr) ) = reshape(u(((idx)*nnf+1):(idx+1)*nnf),fr,fr) ; 
            end
          end
          imagesc(xrange,yrange,I); colorbar; axis equal; hold on;
          h4 = mesh(xg,yg,zeros(s+1));
          h4.FaceColor = 'none';
          h4.EdgeColor = 'k';  
          title('y=0'); hold off
          %offset = offset + nelems(3)*nelems(1)*(nelems(2)/2) ;
          subplot(3,3,5);
          for k=1:nelems(3)
            for j=1:nelems(2)
                idx = offset + (k-1)*(nelems(2)+1)*nelems(1) + (nelems(2)/2)*(nelems(1))  + j-1;
                I( ((j-1)*fr+1):(j*fr),((k-1)*fr+1):(k*fr) ) = reshape(u(((idx)*nnf+1):(idx+1)*nnf),fr,fr) ; 
            end
          end
          imagesc(xrange,yrange,I); colorbar; axis equal; hold on;
          h5 = mesh(xg,yg,zeros(s+1));
          h5.FaceColor = 'none';
          h5.EdgeColor = 'k';  
          title('y=mid'); hold off
          %offset = offset + nelems(3)*nelems(1)*(nelems(2)/2) ;
          subplot(3,3,6);
          for k=1:nelems(3)
            for j=1:nelems(2)
                idx = offset + (k-1)*(nelems(2)+1)*nelems(1) + (nelems(2))*(nelems(1))  + j-1;
                I( ((j-1)*fr+1):(j*fr),((k-1)*fr+1):(k*fr) ) = reshape(u(((idx)*nnf+1):(idx+1)*nnf),fr,fr) ; 
            end
          end
          imagesc(xrange, yrange,I); colorbar; axis equal;hold on;
          h6 = mesh(xg,yg,zeros(s+1));
          h6.FaceColor = 'none';
          h6.EdgeColor = 'k';  
          title('y=end'); hold off
          offset = offset + nelems(3)*nelems(1)*(nelems(2)+1) ;
          %% z=const faces
          I = zeros(nelems(1)*fr, nelems(2)*fr); 
          % z = 0
          subplot(3,3,7);
          for k=1:nelems(2)
            for j=1:nelems(1)
                idx = offset + (k-1)*nelems(1) + j-1;
                I( ((j-1)*fr+1):(j*fr),((k-1)*fr+1):(k*fr) ) = reshape(u(((idx)*nnf+1):(idx+1)*nnf),fr,fr) ;  
            end
          end
          imagesc(xrange,yrange,I); colorbar; axis equal;hold on;
          h7 = mesh(xg,yg,zeros(s+1));
          h7.FaceColor = 'none';
          h7.EdgeColor = 'k';  
          title('z=0'); hold off
          %offset = offset + nelems(1)*nelems(2)*(nelems(3)/2) ;
          subplot(3,3,8);
          for k=1:nelems(2)
            for j=1:nelems(1)
                idx = offset + (k-1)*nelems(1) + j-1 + nelems(2)*nelems(1)*nelems(3)/2;
                I( ((j-1)*fr+1):(j*fr),((k-1)*fr+1):(k*fr) ) = reshape(u(((idx)*nnf+1):(idx+1)*nnf),fr,fr) ;  
            end
          end
          imagesc(xrange,yrange,I); colorbar; axis equal; hold on;
          h8 = mesh(xg,yg,zeros(s+1));
          h8.FaceColor = 'none';
          h8.EdgeColor = 'k';  
          title('z=mid'); hold off
          %offset = offset + nelems(1)*nelems(2)*(nelems(3)/2) ;
          subplot(3,3,9);
          for k=1:nelems(2)
            for j=1:nelems(1)
                idx = offset + (k-1)*nelems(1) + j-1 + nelems(2)*nelems(1)*nelems(3);
                I( ((j-1)*fr+1):(j*fr),((k-1)*fr+1):(k*fr) ) = reshape(u(((idx)*nnf+1):(idx+1)*nnf),fr,fr) ; 
            end
          end
          imagesc(xrange,yrange,I); colorbar; axis equal; hold on;
          h9 = mesh(xg,yg,zeros(s+1));
          h9.FaceColor = 'none';
          h9.EdgeColor = 'k';  
          title('z=end')
        end

        function vol = solve_leaf(grid, u)
            num_elem = grid.Mesh.nelems;
            nnf = grid.Mesh.refel.nnf;
            nr = grid.Mesh.refel.p;

            for k=1:num_elem(3)
                for j=1:num_elem(2)
                    for i=1:num_elem(1)
                        e = sub2ind (num_elem, i, j, k);
                        fid_p = grid.Mesh.get_global_faces(e);
                        % use solution operator to get interior faces
                        bdy_idx = [];
                        % copy trace to element boundaries 
                        for f=1:length(fid_p)
                           bdy_idx = [bdy_idx, ((fid_p(f)-1)*nnf+1):(fid_p(f)*nnf)];
                        end
                        uv = grid.Mesh.S{e}*[u(bdy_idx); 1];

                        vol(((i-1)*(nr-1)+1):i*(nr-1)+1, ((j-1)*(nr-1)+1):j*(nr-1)+1, ((k-1)*(nr-1)+1):k*(nr-1)+1) = reshape(uv, [nr, nr, nr]);
                    end %i
                end %j 
            end %k
        end % solve leaf 

        function r = fas(grid, rc)
            num_elem = grid.Mesh.nelems;
            nnf = grid.Mesh.refel.nnf;
            forder = [1 3 2 5 4 6];
            % for all elems 
            r = grid.get_u0();
            for k=1:num_elem_c(3)
                for j=1:num_elem_c(2)
                    for i=1:num_elem_c(1)
                        ep = sub2ind (num_elem_c, i, j, k);
                        fid_p = grid.Coarse.Mesh.get_global_faces(ep);
                        % use solution operator to get interior faces
                        bdy_idx = [];
                        % copy trace to element boundaries 
                        for f=1:length(fid_p)
                           bdy_idx = [bdy_idx, ((fid_p(f)-1)*nnf+1):(fid_p(f)*nnf)];
                        end
                        up = grid.Coarse.Mesh.S{ep}*[rc(bdy_idx); 1];
                        upf = grid.Mesh.refel.Pint * up;
                        %------------------------------ 
                        ce(1) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k-1);
                        ce(2) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k-1);
                        ce(3) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k-1);
                        ce(4) = sub2ind (num_elem_f, 2*i, 2*j, 2*k-1);
                        ce(5) = sub2ind (num_elem_f, 2*i-1, 2*j-1, 2*k);
                        ce(6) = sub2ind (num_elem_f, 2*i, 2*j-1, 2*k);
                        ce(7) = sub2ind (num_elem_f, 2*i-1, 2*j, 2*k);
                        ce(8) = sub2ind (num_elem_f, 2*i, 2*j, 2*k);
                        for c=1:8
                            fid_c = grid.Mesh.get_global_faces(ce(c));
                            for f=1:3 
                                % 3 shared faces with parent
                                ff = forder(grid.mergemap.ploc(c,f));
                                idx = grid.Mesh.refel.pcf_idx{c, grid.mergemap.ploc(c,f)};
                                r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf)) = r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf)) + upf(idx);
                                % 3 shared faces with siblings
                                ff = forder(grid.mergemap.sloc(c,f));
                                idx = grid.Mesh.refel.pcf_idx{c,grid.mergemap.sloc(c,f)};
                                r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf)) = r(((fid_c(ff)-1)*nnf+1):(fid_c(ff)*nnf)) + upf(idx); 
                            end
                        end
                        
                    end % i
                end % j
            end % k
            % boundaries ?
        end % function fas 

    end %methods
    
end %classdef


% nnf = grid.Mesh.refel.nnf;
%               fr = sqrt(nnf);
%               I = zeros(32,32);
%               for k=1:8
%                 for j=1:8
%                   I( ((k-1)*fr+1):(k*fr), ((j-1)*fr+1):(j*fr)) = reshape(res( ((8*(k-1) + j-1)*nnf+1):(8*(k-1) + j)*nnf),fr,fr) ; 
%                 end
%               end
%               imagesc(I); colorbar; axis equal;
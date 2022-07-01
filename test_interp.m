% test 1d and 2d interpolation operators ...

order = 15;
r = mgps.refel(3, order);

%===========================

sx = sin(0.5*pi.*r.r);
sy = sin(pi.*r.r);

% plot(r.r,sx,'-','LineWidth', 4);hold on
% set(gcf,'unit','norm','position',[0.2 0.2 0.8 0.8])

% parent 

sxp = sx(2:order)
syp = sy(2:order)
% plot(r.r(2:order), sxp, 'ro','MarkerSize', 25,'MarkerFaceColor','red');

% child 
% sxc = r.p_h_1d*sxp;
% plot(r.r_hby2, sxc, '-gs','MarkerSize', 25,'MarkerFaceColor','green');

%%---- test 

% for i=1:(r.p-2)
%     Vri(i,:)     = mgps.basis.polynomial (r.r(2:order), 0, 0, i-1);
%     Vr(i,:)      = mgps.basis.polynomial (r.r, 0, 0, i-1);
% end
% 
% P      = transpose (Vri \ Vr); %flipud(p_h_1d);
% sxc = P*sxp;
% plot(r.r, sxc, '-gs','MarkerSize', 25,'MarkerFaceColor','green');

% 2d 
sx2 = kron(sxp, syp);
figure,subplot(2,3,1);
set(gcf,'unit','norm','position',[0.2 0.2 0.8 0.8])
imagesc(reshape(sx2,order-1,order-1),[-1,1]); axis equal

%% --- prolongation ---

% p1 = r.p_h_1d(1:end/2,:);
% p2 = r.p_h_1d(end/2+1:end,:);
% Ph{1} = kron(p1, p1);
% Ph{3} = kron(p2, p1);
% Ph{2} = kron(p1, p2);
% Ph{4} = kron(p2, p2);
% 
% P = kron(r.p_h_1d, r.p_h_1d);
% sxc2 = P*sx2;
% subplot(2,3,2);
% imagesc(reshape(sxc2,2*(order-1),2*(order-1)),[-1,1]); axis equal

sxc2 = r.Pp * sx2;
subplot(2,3,2);
imagesc(reshape(sxc2,2*(order-1),2*(order-1)),[-1,1]); axis equal
colorbar




%% --- restriction ---

% r1 = r.r_h_1d(:,1:end/2);
% r2 = r.r_h_1d(:,end/2+1:end);
% Rh{1} = kron(r1, r1);
% Rh{3} = kron(r2, r1);
% Rh{2} = kron(r1, r2);
% Rh{4} = kron(r2, r2);
% 
% R = kron(r.r_h_1d, r.r_h_1d);
% sxp2 = R*sxc2;
% subplot(2,3,3);
% imagesc(reshape(sxp2,(order-1),(order-1))); axis equal

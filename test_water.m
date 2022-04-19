k=4;
rhs = @(x,y,z) sin(k*pi*x).*sin(k*pi*y).*sin(k*pi*z);
test_op;
num_elems = [32,32,32];
order = 7;

nlin = num_elems.*order +1;
mu = ones(nlin).*4;

mb = nlin./5;
me = mb.*4;

mu(mb(1):me(1), mb(2):me(2), mb(3):me(3)) = 80;

% smooth mu 
H = fspecial3('ellipsoid',[7,7,7]);
muSmooth = imfilter(mu, H,'replicate');

% sview(muSmooth)
imagesc(muSmooth(:,:,100)); axis equal; colorbar;

